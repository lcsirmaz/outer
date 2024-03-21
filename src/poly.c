/** poly.c  --  combinatorial part using double description method **/

/***********************************************************************
 * This code is part of OUTER, a linear multiobjective problem solver.
 *
 * Copyright (C) 2024 Laszlo Csirmaz, https://github.com/lcsirmaz/outer
 *
 * This program is free, open-source software. You may redistribute it
 * and/or modify under the terms of the GNU General Public License (GPL).
 *
 * There is ABSOLUTELY NO WARRANTY, use at your own risk.
 ***********************************************************************/

/*
* This part of OUTER is the combinatorial vertex enumeration algorithm.
* We maintain an outer approximation of the final polytope. In each stage
* a new vertex of the approximation is chosen, and the oracle is asked
* whether it is part of the final polytope. If yes, it is marked as
* "final". If not, the oracle returns a facet of the final polytope
* separating the point and the polytope, and cuts the approximation.
*
* The algorithm starts with a d dimensional simplex which has d ideal
* vertices at the positive end of the axes, and a "real" vertex. Each
* facet (final and intermediate) has non-negative normal.
*
* Both the vertices and the facets of the outer approximation are stored
* together with adjacency lists which are handled as bitmaps. Using
* bitmaps allows fast operation on vertices and facets.
*/         

#include <stdio.h>
#include <stdint.h>	/* uint32_t, uint64_t */
#include <stdlib.h>
#include <string.h>
#include "main.h"
#include "report.h"
#include "poly.h"
#include "params.h"

#ifdef USETHREADS
/************************************************************************
* Threads
*
* Creating new vertices is split among PARAMS(Threads) threads. The 
* thread logic is the following:
*   MAIN THREAD:
*        loop {
*            ...
*            distribute work (or ask threads to terminate)
*            -- ThreadBarrierForking--
*            do 1/n of the parallel work
*            -- ThreadBarrierJoining --
*            parallel work is done
*            ...
*        }
*  EXTRA THREADS:
*        loop {
*            ...
*            -- ThreadBarrierForking --
*            optionally exit if instructed
*            get the work and do it
*            -- ThreadBarrierJoining --
*        }
*/

#include <pthread.h>
#include <sys/sysinfo.h>  /* get_nprocs() */

/************************************************************************
*
*       T H R E A D S
*
*************************************************************************
*
* Variables used by threads
*
* thread_data_t ThreadData[MAX_THREADS]
*    data block containing data on a single thread */
typedef struct {
    int		id;	/* thread number */
    int		quit;	/* set to 1 if the thread should stop */
    pthread_t	obj;	/* the thread itself */
} thread_data_t;

static thread_data_t ThreadData[MAX_THREADS]; // 0 is unused

/* thread job, the argument is zero for the calling thread */
typedef void ThreadJob_t(int threadId);
static ThreadJob_t *ThreadJob;
#define ThreadNo	PARAMS(Threads)

/* pthread_barrier_t ThreadBarrierForking, ThreadBarrierJoining
*    barriers for synchronizing work */
static pthread_barrier_t ThreadBarrierForking;
static pthread_barrier_t ThreadBarrierJoining;

/* void *extra_thread(void *arg)
*    the code for extra threads - forward declaration */
static void *extra_thread(void *arg);

/* int create_threads(void)
*    create barriers, mutex, threads, and start them. Threads will
*    wait until they get their first assignment. */
int create_threads(void){
int i,rc;
    if(PARAMS(ProblemObjects)<2){
        PARAMS(Threads)=1;
        return 0; // do not use threads in this case
    }
    if(PARAMS(Threads)<=0){ // figure out the number of CPU's
        PARAMS(Threads)=get_nprocs();
    }
    if(PARAMS(Threads)<=1){
        PARAMS(Threads)=1;
        return 0; // nothing to do
    }
    if(PARAMS(Threads)>MAX_THREADS){ PARAMS(Threads)=MAX_THREADS; }
    // barriers
    if((rc = pthread_barrier_init(&ThreadBarrierJoining, NULL, PARAMS(Threads)))){
        report(R_fatal,"error: pthread_barrier_init, rc: %d\n", rc);
        return 1;
    }
    if((rc = pthread_barrier_init(&ThreadBarrierForking, NULL, PARAMS(Threads)))){
        report(R_fatal,"error: pthread_barrier_init, rc: %d\n", rc);
        return 1;
    }
    // create and start threads, they will stop arriving at the barriers
    report(R_info,"Using %d threads\n", PARAMS(Threads));
    ThreadData[0].id=0;
    for(i=1;i<PARAMS(Threads);i++){
        ThreadData[i].id = i;
        ThreadData[i].quit = 0;
        if((rc = pthread_create(&(ThreadData[i].obj),NULL,extra_thread,&ThreadData[i]))){
            report(R_fatal,"error creating thread %d, rc: %d\n",i,rc);
            return 1;
        }
    }
    return 0;
}
/* void thread_execute(ThreadJob_t job)
*    execute "job" using all avaiable threads */
static void thread_execute(ThreadJob_t job)
{   if(PARAMS(Threads)<2){ // single thread
        job(0); return;
    }
    ThreadJob=job;
    pthread_barrier_wait(&ThreadBarrierForking); // start threads
      job(0); // main thread
    pthread_barrier_wait(&ThreadBarrierJoining); // wait until others finish
}

/* void *extra_thread(void *arg)
*    each thread waits, then calls the routine in ThreadJob, then sleeps */
static void *extra_thread(void *arg){
thread_data_t *data = (thread_data_t*)arg;
int myId=data->id;
    report(R_info,"Thread %d started ...\n",myId);
    while(1){
       pthread_barrier_wait(&ThreadBarrierForking);
       if(data->quit) break; // stop
       ThreadJob(myId);
       pthread_barrier_wait(&ThreadBarrierJoining);
    }
    report(R_info,"Thread %d stopped\n",myId);
    return NULL;
}
/* void stop_threads(void)
*    tell extra threads to exit and wait for (reap) them.
*    This blocks until all threads exit. */
void stop_threads(void){
int i;
    if(PARAMS(Threads)<2) return;
    report(R_info,"Stopping threads ...\n");
    for(i=1;i<PARAMS(Threads);i++){
        ThreadData[i].quit = 1;
    }
    pthread_barrier_wait(&ThreadBarrierForking);
    for(i=1;i<PARAMS(Threads);i++){
        pthread_join(ThreadData[i].obj,NULL);
    }
}
#else /* ! USETHREADS */
#define ThreadNo		1 /* number of threads */
#endif /* USETHREADS */

/************************************************************************
*
*    M E M O R Y   M A N A G E M E N T 
*
*************************************************************************
*
* Memory is arranged in slots; a slot consists of several blocks of
* given blocksize. When reallocating, either the number of blocks, the
* blocksize, or both can change. When blocksize changes, the old content
* is adjusted to match the new blocksize.
* Content in the main slots are needed for the duration of the algorithm,
* and the blocksize can change. Temporary slots are used within a single
* iteration, and blocksize is fixed.
*
* void report_memory_usage(channel,prompt)
*   report the cnages since the last report to the given channel
*
* type *get_memory_ptr(type,slot)
*   retrieve the first memory block in the given slot
*
* void clear_memory_slote(void)
*   initialize all memory slots, should be called first
*
* void yalloc(type,slot,blocks,blocksize)
*   initialize a main slot; each block has blocksize words of type
*
* void yrequest(type,slot,blocks,blocksize)
*   request memory allocation change in the given slot. Should be
*   followed by calling reallocmem()
*
* int reallocmem(void)
*   perform previously requested memory allocation changes.
*
* void yfree(slot)
*   release all memory from the given slot
*
* void talloc(type,slot,blocks,blocksize)
*   initialize a temporary slot; previous content is released first.
*   use only for global memory allocation
*
* void talloc2(type,slot,threadID,blocks,blocksize)
*   same as talloc for private memory allocation
*
* void trequest(slot,threadID,blocks)
*   reallocate private memory to 'blocks' many blocks
*/

typedef enum {	/* main memory slots */
M_VertexCoordStore	= 0,	/* vertex coordinates */
M_FacetCoordStore,		/* facet coordinates */
M_VertexAdjStore,		/* adjacency list of vertices */
M_FacetAdjStore,		/* adjacency list of facets */
M_VertexLiving,			/* single vertex bitmap of actual vertices */
M_VertexFinal,			/* single vertex bitmap of final vertices, subset of VertexLiving */
M_MAINSLOTS,			/* last main slot index */
		/* temporary slots - global for all threads */
M_VertexDistStore=M_MAINSLOTS,	/* vertex distances from the new facet */
M_VertexPosnegList,		/* indices of vertices on positive/negative side */
		/* private slots for threads */
M_THREAD_SLOTS,
M_FacetList=M_THREAD_SLOTS,	/* facets adjacent to two vertices */
M_VertexWork,			/* facet bitmap */
M_FacetArray,			/* calculating vertex coordinates */
M_NewVertexCoordStore,		/* new vertex coordinates */
M_NewVertexAdjStore,		/* adjacency list of new vertices */

M_THREAD_RESERVED_PLACE,	/* memory block for the next thread */
M_THREAD_RESERVED_PLACE_END = M_THREAD_SLOTS+((M_THREAD_RESERVED_PLACE-M_THREAD_SLOTS)*MAX_THREADS)-1,

M_MSLOTSTOTAL			/* total number of managed memory */
} memslot_t;

/* int NUM_M_THREAD_SLOTS
*    number of thread-specific memory slots */
#define NUM_M_THREAD_SLOTS	(M_THREAD_RESERVED_PLACE-M_THREAD_SLOTS)

/* memslot_t M_thread(slot,threadId)
*    thread specific slot id assigned to the given thread */
#define M_thread(slot,threadId)		\
    ((slot)+(threadId)*NUM_M_THREAD_SLOTS)

/* slot title for reporting */
#define TM_VertexCoordStore	"VertexCoord"
#define TM_FacetCoordStore	"FacetCoord"
#define TM_VertexAdjStore	"VertexAdj"
#define TM_FacetAdjStore	"FacetAdj"
#define TM_VertexLiving		"VertexLiving"
#define TM_VertexFinal		"VertexFinal"
#define TM_VertexDistStore	"VertexDist"
#define TM_VertexPosnegList	"VertexPosNeg"
#define TM_FacetList		"FacetList"
#define TM_VertexWork		"VertexWork"
#define TM_FacetArray		"FacetArray"
#define TM_NewVertexCoordStore	"NewVertexCoord"
#define TM_NewVertexAdjStore	"NewVertexAdj"

/* MEMSLOT
*    the memory slot structure */

typedef struct { /* memory slot structure */
  size_t blocksize;	/* block size (in bytes) */
  size_t blockno;	/* number of blocks */
  size_t newblocksize;	/* new block size */
  size_t newblockno;	/* new block count */
  size_t rsize;		/* real size in bytes */
  void   *ptr;		/* pointer to block zero */
  const char *title;	/* block title */
  const char *type;	/* basic type */
  size_t bsize;		/* size of item type */
  size_t rreport;	/* rsize at last report */
} MEMSLOT;

/* struct MEMSLOT memory_slots[]
*    static array containing for each slot the actual blocksize, number
*    of blocks, and a pointer to the actual location. The location can
*    change when reallocating any other memory slot.
*/
static MEMSLOT memory_slots[M_MSLOTSTOTAL]; /* memory slots */

/* bool OUT_OF_MEMORY
*    flag indicating whether we are out of memory.
*/
#define OUT_OF_MEMORY	dd_stats.out_of_memory

/* type *get_memory_ptr(type,slot)
*    the actual memory block in the given slot.
*      type:  storage type of the items in the block
*      slot:  the memory slot
*/
#define get_memory_ptr(type,slot)	\
    ((type *)memory_slots[slot].ptr)

/* void report_memory_usage(channel,prompt)
*    for each used slot report the blocksize, number of blocks,
*    total memory used by this slot, and the actual pointer.
*      channel: report type channel
*      prompt:  prompt for the first line
*/
void report_memory_usage(report_type ch,int force,const char *prompt)
{int i; MEMSLOT *ms; char buff[50];
    for(i=0,ms=&memory_slots[0];i<M_MSLOTSTOTAL;i++,ms++) 
    if(ms->ptr && ms->bsize>0 && (force!=0 || ms->rsize != ms->rreport) ) {
        ms->rreport=ms->rsize;
        if(prompt){
            report(ch,"%s\n          ---slot----------+--size---+---blocks--+---blocksize---\n",prompt);
            prompt=NULL;
        }
        if(ms->rsize<1000ul){ sprintf(buff,"%zu",ms->rsize); }
        else if(ms->rsize<1000000ul){ sprintf(buff,"%.2fk",(double)(ms->rsize)*1e-3); }
        else if(ms->rsize<1000000000ul){ sprintf(buff,"%.2fM",(double)(ms->rsize)*1e-6); }
        else { sprintf(buff,"%.2fG",(double)(ms->rsize)*1e-9); }
        report(ch,"          %2d %-14s|%8s |%10zu | %zu * %s\n",
           i+1,ms->title,buff,ms->blockno,ms->blocksize/ms->bsize,ms->type);
    }
    if(prompt==NULL){
       report(ch,"          -----------------+---------+-----------+---------------\n");
    }
}

/* void clear_memory_slots(void)
*    clear all entries in all slots */
#define clear_memory_slots()	\
    memset(&memory_slots[0],0,sizeof(memory_slots))


/* void yalloc(type,slot,n,bsize)
*    initializes a main memory slot by requesting n blocks; each
*    block is an array of bsize elements of the given type.
       type:  storage type of an item
       slot:  memory slot
       n:     number of blocks
       bsize: number of items in a block
*/
#define yalloc(type,slot,n,bsize)	\
    AUX_init_main_slot(slot,n,bsize,sizeof(type),T##slot,mkstringof(type))

/* void AUX_init_main_slot(slot,nno,n,bsize,title,type)
*    initializes a main slot by allocating the requested memory.
*     slot:    memory slot to be initialized
*     nno:     number of initial blocks
*     n:       number of items in a block
*     bsize:   size of an item (in bytes)
*     title:   memory slot title (string) for reporting
*     type:    item type (string) for reporting
*/
static void AUX_init_main_slot(memslot_t slot, size_t nno, size_t n, size_t bsize,
                           const char *title, const char *type)
{size_t total,nsize; MEMSLOT *ms;
    if(OUT_OF_MEMORY) return;
    ms=&memory_slots[slot];
    ms->newblocksize=0;
    ms->newblockno=0;
    ms->rreport=0;
    ms->bsize=bsize;	// item size in bytes
    ms->title=title;	// slot title
    ms->type=type;	// item type as string
    nsize=n*bsize;	// block size in bytes
    total=nno*nsize;	// total requested size in bytes
    ms->blocksize=nsize;
    ms->blockno=nno;	// number of requested blocks
    ms->rsize=total;
    ms->ptr=malloc(total);
    if(!ms->ptr){
        report(R_fatal,
           "Out of memory for slot=%d (%s), blocksize=%zu, n=%zu\n",
           slot,title,nsize,nno);
        OUT_OF_MEMORY=1;
        return;
    }
    dd_stats.total_memory += total;
    if(dd_stats.total_memory>dd_stats.max_memory)
            dd_stats.max_memory=dd_stats.total_memory;
    memset(ms->ptr,0,total);
    return;
}

/* void yrequest(type,slot,n,bsize)
*    request memory at a main slot; should be followed by calling
*    reallocmem().
*       type:  storage type of an item
*       slot:  memory slot
*       n:     number of blocks requested
*       bsize: number of items in a block
*/
#define yrequest(type,slot,n,bsize)	\
    AUX_request_main_mem(slot,n,(bsize)*sizeof(type))

/* void AUX_request_main_mem(slot,nno,nsize)
*    records the requested block count and block size. The actual 
*    memory allocation is done by reallocmem()
*      slot:  memory slot
*      nno:   number of blocks requested
*      nsize: size of a block (in bytes)
*/
static inline void AUX_request_main_mem(memslot_t slot, size_t nno, size_t nsize)
{MEMSLOT *ms;
    ms=&memory_slots[slot];
    if(ms->blocksize==nsize && nno<=ms->blockno) return;
    ms->newblocksize=nsize;
    ms->newblockno=nno;
}    

/* int reallocmem(void)
*    allocate previously requested memory; if successful, adjust
*    blocks to the given count and size, and clear the new part.
*    Return 1 if out of memory, otherwise return 0. */
static int reallocmem(void)
{MEMSLOT *ms; int j,success; size_t total; void *ptr;
    for(j=0,success=1,ms=&memory_slots[0];
        success && j<M_MAINSLOTS; j++,ms++
    ){
        if(ms->newblocksize){
            total=ms->newblocksize*ms->newblockno;
            if(ms->rsize<total){
                dd_stats.memory_allocated_no++;
                ptr=realloc(ms->ptr,total);
                if(ptr){
                    ms->ptr=ptr; 
                    dd_stats.total_memory += total-ms->rsize;
                    if(dd_stats.total_memory>dd_stats.max_memory)
                          dd_stats.max_memory=dd_stats.total_memory;
                    ms->rsize=total;
                }
                else { success=0; }
            }
        }
    }
    if(!success){ OUT_OF_MEMORY=1; return 1; }
    // adjust block structure
    for(j=0,ms=&memory_slots[0]; j<M_MAINSLOTS; j++,ms++){
        if(ms->newblocksize){
            if(ms->blocksize < ms->newblocksize){ // block size grow
                size_t i,n; char *bo,*bn; // pointer to old and new
                n=ms->newblockno; if(ms->blockno<n) n=ms->blockno;
                bo=((char*)ms->ptr)+(n*ms->blocksize);
                bn=((char*)ms->ptr)+(n*ms->newblocksize);
                for(i=n;i>0;i--){
                    bo -= ms->blocksize; bn -= ms->newblocksize;
                    if(i>1) memmove(bn,bo,ms->blocksize);
                    // clear the rest
                    memset(bn+ms->blocksize,0,ms->newblocksize-ms->blocksize);
                }
            } else if(ms->blocksize > ms->newblocksize) { // block size shrank
                size_t i,n; char *bo,*bn; // pointers to old and new
                n=ms->newblockno; if(ms->blockno<n)n=ms->blockno;
                bo=(char*)ms->ptr; bn=(char*)ms->ptr;
                for(i=0;i<n;i++){
                    if(i) memmove(bn,bo,ms->newblocksize);
                    bo += ms->blocksize;
                    bn += ms->newblocksize;
                }
            }
            ms->blocksize=ms->newblocksize;
            if(ms->blockno < ms->newblockno){ // clear the rest
                memset(((char*)ms->ptr)+ms->blockno*ms->blocksize,0,
                     (ms->newblockno-ms->blockno)*ms->blocksize);
            }
            ms->blockno=ms->newblockno;
            ms->newblocksize=0;
            ms->newblockno=0;
    // shrink too large allocations
            total=ms->blocksize*ms->blockno;
            if(ms->rsize>total+DD_HIGHWATER){
                ptr=realloc(ms->ptr,total+DD_LOWWATER);
                if(ptr){
                    ms->ptr=ptr;
                    dd_stats.total_memory -= ms->rsize-(total+DD_LOWWATER);
                    ms->rsize=total+DD_LOWWATER;
                }
            }
        }
    }
    return 0;
}

/* void yfree(slot)
*    free the memory in the given slot */
inline static void yfree(memslot_t slot)
{MEMSLOT *ms;
    ms=&memory_slots[slot];
    ms->newblocksize=0;
    ms->newblockno=0;
    ms->blockno=0;
    ms->rsize=0;
    if(ms->ptr) free(ms->ptr);
    ms->ptr=(void*)0;
}

/* void talloc(type,slot,n,bsize)
*    request initial memory for a temporary slot. The allocated memory
*    is not cleared
*        type:  storage type of an item
*        slot:  memory slop
*        n:     number of requested blockss
*        bsize: size of an item in bytes
*  void talloc2(type,slot,thread_ID,n,bsize)
*    talloc() for private slot for a thread */
#define talloc(type,slot,n,bsize)	talloc2(type,slot,0,n,bsize)
#define talloc2(type,slot,threadID,n,bsize)	\
    AUX_init_temp_slot(M_thread(slot,threadID),n,bsize,\
           sizeof(type),T##slot,mkstringof(type))

/* void AUX_init_temp_slot(slot,nno,n,bsize,title,type)
*    requests nno blocks, each of size n*bsize at the given slot.
*    The memory is not cleared; should check OUT_OF_MEMORY
*      slot:    memory slot to be initialized
*      nno:     number of requested blocks
*      n:       number of items in each block
*      bsize:   size of an item in bytes
*      title:   memory slot title (string) for reporting
*      type:    item type (string) for reporting 
*/
static void AUX_init_temp_slot(memslot_t slot, size_t nno, size_t n, size_t bsize,
                           const char *title, const char *type)
{size_t total,nsize; MEMSLOT *ms;
    if(OUT_OF_MEMORY) return;
    ms=&memory_slots[slot];
    nsize=n*bsize;
    if(ms->bsize!=bsize || ms->blocksize!=nsize || ms->blockno<nno){
        ms->bsize=bsize; ms->blocksize=nsize; ms->blockno=nno;
    }
    ms->title=title; ms->type=type;
    total=nno*nsize;
    if(total <= ms->rsize) return;
    if(ms->ptr){ 
        free(ms->ptr);
        dd_stats.total_memory -= ms->rsize;
        dd_stats.memory_allocated_no++;
    }
    ms->rsize=total;
    dd_stats.total_memory += total;
    if(dd_stats.total_memory>dd_stats.max_memory)
          dd_stats.max_memory=dd_stats.total_memory;
    ms->ptr = malloc(total);
    if(!ms->ptr){
        report(R_fatal,"Out of memory for slot=%d (%s), blocksize=%zu, n=%zu\n",
            slot,title,nsize,nno);
        OUT_OF_MEMORY=1;
    }
}

/* void trequest(slot,threadID,n)
*    expand the temporary slot to n blocks, keep the previous block
*    size. The new memory is not cleared. */
#define trequest(slot,threadID,n)	\
    AUX_request_temp_mem(M_thread(slot,threadID),n)

/* void AUX_request_temp_mem(slot,nno)
*    request more blocks for the initialized temporary memory slot
*       slot:  memory slot
*       nno:   number of blocks required
*/
static inline void AUX_request_temp_mem(memslot_t slot,size_t nno)
{size_t total; MEMSLOT *ms; void *ptr;
    if(OUT_OF_MEMORY) return;
    ms=&memory_slots[slot];
    total = ms->blocksize*nno;
    if(total<=ms->rsize) return;
    dd_stats.memory_allocated_no++;
    ptr=realloc(ms->ptr,total);
    if(!ptr){ OUT_OF_MEMORY=1; return; }
    dd_stats.total_memory += total - ms->rsize;
    if(dd_stats.total_memory>dd_stats.max_memory)
          dd_stats.max_memory=dd_stats.total_memory;
    ms->rsize=total;
    ms->blockno=nno;
    ms->ptr=ptr;
}

/* void free_adjacency_lists(void)
*    after a break request, vertex and facet adjancy lists are not used
*    release their memory */
void free_adjacency_lists(void)
{   yfree(M_VertexAdjStore); yfree(M_FacetAdjStore); }

/************************************************************************
*
*     B I T M A P S
*
*************************************************************************
*
* A *bitmap* is an array of 64 bit (or 32 bit) unsigned integers
* accommodating the requested number of bits. Bitmaps are used as 
* adjacency lists, and also to store vertex and facet attributes.
*
* BITMAP_t, BITMAP0, BITMAP1, packshift, packmask
*   Bitmaps are declared as type BITMAP_t *B. BITMAP0 and BITMAP1
*   are zero and 1 in BITMAP_t type. To get/set "index" bit in
*   bitmap B use packshift and packmask as follows:
*     (B[index>>packshift]>>(index&packmask)) & 1
*     B[index>>packshift] |= BITMAP1 << (index&packmask)
*
* bool extract_bit(from,cnt)
*    extract bit cnt from a BITMAP_t array
*
* void clear_bit(from,cnt)
*    clear bit cnt from a BITMAP_t array
*
* void set_bit(from,cnt)
*    set bit cnt in the BITMAP_t array
*
* int get_bitcount(BITMAP_t v)
*    count the bits set in a single BITMAP_t word
*/

#ifdef BITMAP_32		/* 32 bit bitmap */
typedef uint32_t BITMAP_t;
#define packsizelog	2	/* sizeof(BITMAP_t)== 1<<packsizelog */
#else				/* 64 bit bitmap */
typedef uint64_t BITMAP_t;
#define packsizelog	3	/* sizeof(BITMAP_t)== 1<<packsizelog */
#endif /* BITMAP_32 */

#define packsize	(1<<packsizelog)	/* 4 or 8 bytes */
#define packshift	(packsizelog+3) 	/* in bits (5 or 6 )*/
#define packmask	((1<<packshift)-1)	/* 31 or 63 */
#define BITMAP1		((BITMAP_t)1u)    	/* BITMAP_t unsigned 1 */
#define BITMAP0		((BITMAP_t)0u)		/* BITMAP_t unsigned zero */

/* bool extract_bit(BITMAP_t bm[], int cnt)
*  void clear_bit  (BITMAP_t bm[], int cnt)
*  void set_bit    (BITMAP_t bm[], int cnt) */
#define extract_bit(bm,cnt)	\
    (((bm)[(cnt)>>packshift]>>((cnt)&packmask))&1)

#define clear_bit(bm,cnt)	\
    (bm)[(cnt)>>packshift] &= ~(BITMAP1<<((cnt)&packmask))

#define set_bit(bm,cnt)		\
    (bm)[(cnt)>>packshift] |= BITMAP1<<((cnt)&packmask)

/* int get_bitcount(BITMAP_t v)
*    cont the number of bits in v */
#ifndef NO_ASM
/* this routine uses the population count (popcnt) machine code
   to get the number of bits in a word. */
static inline int get_bitcount(BITMAP_t v)
{register BITMAP_t res;
    asm ("popcnt %[w], %[t]"
         :[t] "=rm" (res)
         :[w] "rm"  (v));
    return (int)res;
}

#else /* no assembler code */
static int get_bitcount(BITMAP_t v)
{
//    these exceptional cases do not seem to help
//    if(v==BITMAP0){ return 0; }
//    if(v==~BITMAP0){ return 1<<packshift; }
#   ifdef BITMAP_32	/* 32 bit bitmap */
    v=(v&0x55555555u)+((v>>1)&0x55555555u);
    v=(v&0x33333333u)+((v>>2)&0x33333333u);
    v=(v&0x0F0F0F0Fu)+((v>>4)&0x0F0F0F0Fu);
    v=(v&0x00FF00FFu)+((v>>8)&0x00FF00FFu);
    return (int)((v&0x0000FFFF)+(v>>16));
#   else		/* 64 bit bitmap */
    v=(v&0x5555555555555555ul)+((v>>1)&0x5555555555555555ul);
    v=(v&0x3333333333333333ul)+((v>>2)&0x3333333333333333ul);
    v=(v&0x0F0F0F0F0F0F0F0Ful)+((v>>4)&0x0F0F0F0F0F0F0F0Ful);
    v=(v&0x00FF00FF00FF00FFul)+((v>>8)&0x00FF00FF00FF00FFul);
    v=(v&0x0000FFFF0000FFFFul)+((v>>16)&0x0000FFFF0000FFFFul);
    return (int)((v&0xFFFF)+(int)(v>>32));
#   endif
}
#endif /* NO_ASM */

/*************************************************************************
*
*    S T A T I S T I C S
*
**************************************************************************/

DD_STATS dd_stats;

/*************************************************************************
*
*    V E R T I C E S    A N D    F A C E T S
*
**************************************************************************
*
* int DIM
*   problem dimension = number of objectives
* int VertexSize, NextVertex, MaxVertices, VertexBitmapBlockSize
* int FacetSize, NextFacet, MaxFacets, FacetBitmapBlockSize
*   size of a blocks; next free block index; maximal available
*   index; size of the corresponding bitmap block
* int ThisFacet
*   index of facet we are adding to the approximation
* double *VertexCoords(vno), BITMAP_t *VertexAdj(vno)
* double *FacetCoords(fno), BITMAP_t *FacetAdj(fno)
*   the memory block and adjacency block of a vertex and a facet
*
* BITMAP_t *VertexLiving, *VertexFinal
*   bitmaps marking valid and final vertices
*
* double VertexDist[0 .. MaxVertices]
*   distance of vertices from the next facet
*
* int *VertexPosnegList
*   vertices on the positive and negative side of the new facet
*/

#define DIM	PARAMS(ProblemObjects)

static int
  VertexSize,		// size of of a vertex block = DIM
  NextVertex,		// next free index for a vertex
  MaxVertices,		// number of bits in a facet bitmap
  VertexBitmapBlockSize,// size of the vertex bitmap block
  FacetSize,		// size of the facet equation block = DIM+1
  NextFacet,		// next free slot for a facet
  MaxFacets,		// number of bits in a facet bitmap
  FacetBitmapBlockSize, // size of the facet bitmap block
  ThisFacet;		// the facet we are working with

/* double *VertexCoords(vno); BITMAP_t *VertexAdj(vno);
      coordinate (DIM) and adjacency list (FacetBitmapBlockSize) of 
      vertex 'vno'
   double *FacetCoords(fno); BITMAP_t *FacetAdj(fno)
      coordinate (DIM+1) and adjacency list (VertexBitmapBlockSize)
      of facet 'fno'
*/
#define VertexCoords(vno)	\
    (get_memory_ptr(double,M_VertexCoordStore)+((vno)*VertexSize))
#define VertexAdj(vno)		\
    (get_memory_ptr(BITMAP_t,M_VertexAdjStore)+((vno)*FacetBitmapBlockSize))
#define FacetCoords(fno)	\
    (get_memory_ptr(double,M_FacetCoordStore)+((fno)*FacetSize))
#define FacetAdj(fno)		\
    (get_memory_ptr(BITMAP_t,M_FacetAdjStore)+((fno)*VertexBitmapBlockSize))

/* BITMAP_t *VertexLiving, *VertexFinal */
#define VertexLiving		\
    get_memory_ptr(BITMAP_t,M_VertexLiving)
#define VertexFinal		\
    get_memory_ptr(BITMAP_t,M_VertexFinal)
/* double VertexDist(vno); int *VertexPosNegList */
#define VertexDist(vno)		\
    get_memory_ptr(double,M_VertexDistStore)[vno]
#define VertexPosnegList	\
    get_memory_ptr(int,M_VertexPosnegList)
 
/* void get_dd_vertexno(void)
*    compute the number of living and final vertices intto dd_stats
*  int get_vertexnum()
*    number of livint vertices
*  int get facetnum()
*    number of facets */
void get_dd_vertexno(void)
{int i,vno; BITMAP_t vl,vf;
    dd_stats.living_vertex_no=-DIM; dd_stats.final_vertex_no=-DIM;
    for(i=0;i<VertexBitmapBlockSize;i++){
        vl=VertexLiving[i];
        dd_stats.living_vertex_no += get_bitcount(vl);
        vf=VertexFinal[i];
        dd_stats.final_vertex_no += get_bitcount(vf);
        if(~vl & vf){ // consistency checking
            vno=i<<packshift;
            while(!((~vl&vf)&1)){ vno++;vl>>=1;vf>>=1;}
            report(R_err,"Consistency error: final but not living vertex %d\n",vno);
        }
    }
    if(dd_stats.living_vertex_no<0) dd_stats.living_vertex_no=0;
    if(dd_stats.final_vertex_no<0)  dd_stats.final_vertex_no=0;
}
int get_vertexnum(void)
{int i,total; BITMAP_t v;
    total=-DIM;
    for(i=0;i<VertexBitmapBlockSize;i++){
        v=VertexLiving[i]; total += get_bitcount(v);
    }
    return total;
}
int get_facetnum(void)
{   return NextFacet; }

/*************************************************************************
*
*    I T E R A T I O N
*
**************************************************************************
*
* The main iteration adds a new facet to the actual approximation. First,
* for each vertex 'vno' VertexDist(vno) contains the distance of the vertex
* from the facet. Vertex numbers with positive distances are collected at
* the from of VertexPostNextList, vertex numbers with negative distance are
* collected from the end. For each postive / negative pair (vp,vn) the
* following are computed by each thread:
*  FacetList       -- list of facet numbers which are adjacent to both vp,vn
*  WertexWork      -- BITMAP of vertices adjacent to facets in FacetList
* If (vp,vn) is an edge, it intersects the new facet in a new vertex.
*  NewVertexCoords -- coordinates of new vertices
*  NewVertexAdj    -- adjacency lists of new vertices
*  FacetArray      -- equation of all facets adjacent to the new vertex
* Newly created vertices by each thread:
*  NewVertex       -- number of new vertices generated
*  MaxNewVertex    -- maximal space available to the thread
*/

/* When there are no threads, thId=0 */
#define FacetList(thId)			\
    get_memory_ptr(int,M_thread(M_FacetList,thId))
#define VertexWork(thId)		\
    get_memory_ptr(BITMAP_t,M_thread(M_VertexWork,thId))
/* double *NewVertexCoords(thId,vno); BITMAP_t *NewVertexAdj(ThId,vno); */
#define NewVertexCoords(thId,vno)	\
    (get_memory_ptr(double,M_thread(M_NewVertexCoordStore,thId))+\
     ((vno)*VertexSize))
#define NewVertexAdj(thId,vno)		\
    (get_memory_ptr(BITMAP_t,M_thread(M_NewVertexAdjStore,thId))+\
      ((vno)*FacetBitmapBlockSize))
#define FacetArray(thId)		\
    get_memory_ptr(double,M_thread(M_FacetArray,thId))

static int
  NewVertex[MAX_THREADS],    // number of newly created vertices
  MaxNewVertex[MAX_THREADS], // available space
  ErrorNo[MAX_THREADS];	     // consistency errors

/* BOOL is_livingVertex(vno)
*     check if bit 'vno' is set in bitmap VertexLiving
* void set_in_VertexLiving(vno)
*     set bit 'vno' in bitmap VertexLiving
* BOOL is_finalVertex(vno)
*     check if bit 'vno' is set in bitmap VertexFinal
* void set_in_VertexFinal(vno)
*     set bit 'vno' in bitmap VertexFinal */

#define is_livingVertex(vno)	    extract_bit(VertexLiving,vno)
#define set_in_VertexLiving(vno)    set_bit(VertexLiving,vno)
#define is_finalVertex(vno)	    extract_bit(VertexFinal,vno)
#define set_in_VertexFinal(vno)     set_bit(VertexFinal,vno)

/* void mark_vertex_as_final(vno)
*     exported version of set_in_VertexLiving(vno)
*  void intersect_FacetAdj_VertexLiving(fno)
*     clear bits in the adjacency list of 'fno' which are not living
*  void move_vertex_to(coords,adj)bitmap,vno)
*     move vertex with coordinates and adjacency list to the given index
*  void clear_VertexAdj(vno)
*     clear the adjacecny list of vertex 'vno'
*  void clear_FacetAdj(fno)
*     clear the adjacency list of facet 'fno'
* void copy_VertexLiving_to_(where)
*     copy the bitmap VertexLiving to the given address
* void move_NewVertex_th(threadId,vno)
*      move NewVertex[threadId]-th element of NewVertexCoords() and
*      NewVertexAdj to the index 'vno' */

void mark_vertex_as_final(int vno)
{   set_bit(VertexFinal,vno); }

inline static void intersect_FacetAdj_VertexLiving(int fno)
{int i;
    for(i=0;i<VertexBitmapBlockSize;i++)
        FacetAdj(fno)[i] &= VertexLiving[i];
}
inline static void move_vertex_to(double *coords,BITMAP_t *adj,int vno)
{   memcpy(VertexCoords(vno),coords,VertexSize*sizeof(double));
    memcpy(VertexAdj(vno),adj,FacetBitmapBlockSize*sizeof(BITMAP_t)); }

inline static void clear_VertexAdj(int vno)
{   memset(VertexAdj(vno),0,FacetBitmapBlockSize*sizeof(BITMAP_t)); }

inline static void clear_FacetAdj(int fno)
{   memset(FacetAdj(fno),0,VertexBitmapBlockSize*sizeof(BITMAP_t)); }

#define copy_VertexLiving_to(where)	\
    memcpy(where,VertexLiving,VertexBitmapBlockSize*sizeof(BITMAP_t))

inline static void move_NewVertex_th(int thId,int vno)
{   move_vertex_to(NewVertexCoords(thId,NewVertex[thId]),
                   NewVertexAdj(thId,NewVertex[thId]),vno); }

/************************************************************************
* Initialization
*
* int init_dd_structure(int vertexno,int facetno)
*    initialize the DD algorithm structures. Return value:
*   0  initialization OK
*   1  problem dimension is too large or out of memory
* int init_dd(double first_vertex[0:dim-1]
*    start the DD algoritm with the first external vertex. Calls
*    init_dd_structure(0,0); */

int init_dd_structure(int vertexno,int facetno)
{int i,j;
    if(DIM<1 || DIM>MAXIMAL_ALLOWED_DIMENSION){
        report(R_fatal,"Problem dimension %d is out of allowed range (1 .. %d)\n",
           DIM,MAXIMAL_ALLOWED_DIMENSION);
        return 1;
    }
    MaxVertices=vertexno;
    if(MaxVertices<DD_INITIAL_VERTEXNO) MaxVertices=DD_INITIAL_VERTEXNO;
    MaxFacets=facetno;
    if(MaxFacets<DD_INITIAL_FACETNO) MaxFacets=DD_INITIAL_FACETNO;
    while(MaxVertices<DIM) MaxVertices += (DD_VERTEX_ADDBLOCK<<packshift);
    while(MaxFacets<DIM) MaxFacets += (DD_FACET_ADDBLOCK<<packshift);
    FacetSize=DIM+1; VertexSize=DIM;
    // bitmaps
    VertexBitmapBlockSize = (MaxVertices+packmask)>>packshift;
    FacetBitmapBlockSize = (MaxFacets+packmask)>>packshift;
    // clear statistics data
    memset(&dd_stats,0,sizeof(dd_stats));
    dd_stats.vertices_allocated_no=1;
    dd_stats.vertices_allocated=MaxVertices;
    dd_stats.facets_allocated_no=1;
    dd_stats.facets_allocated=MaxFacets;
    dd_stats.facetno=0; // no facets added, the ideal facet is not stored
    // clear memory slots
    clear_memory_slots();
    yalloc(double,M_VertexCoordStore,MaxVertices,VertexSize); // VertexCoordStore
    yalloc(double,M_FacetCoordStore,MaxFacets,FacetSize); // FacetCoordStore
    yalloc(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize); // VertexAdjStore
    yalloc(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize); // FacetAdjStore
    yalloc(BITMAP_t,M_VertexLiving,1,VertexBitmapBlockSize); // VertexLiving
    yalloc(BITMAP_t,M_VertexFinal,1,VertexBitmapBlockSize); // VertexFinal
    if(OUT_OF_MEMORY) return 1;
    dd_stats.memory_allocated_no=1;
    NextVertex=0; NextFacet=0;
    // ideal vertices with index 0 .. DIM-1
    for(i=0;i<DIM;i++){ // i-th ideal vertex
        for(j=0;j<DIM;j++) VertexCoords(NextVertex)[j]= i==j?1.0: 0.0;
        clear_VertexAdj(NextVertex); // the ideal facet is not stored
        set_in_VertexLiving(NextVertex); 
        set_in_VertexFinal(NextVertex);
        NextVertex++;
    }
    return 0;
}

/* int init_dd(double first_vertex[1:DIM-1])
*    set up the first outer approximation using the given vertex */
void init_dd(const double *coords)
{int i,j;
    // add the vertex of the first approximation
    for(j=0;j<DIM;j++)VertexCoords(NextVertex)[j]=coords[j];
    clear_VertexAdj(NextVertex);
    set_in_VertexLiving(NextVertex);
    // add DIM coordinate facets incident to the given point
    for(i=0;i<DIM;i++){
        for(j=0;j<DIM;j++)FacetCoords(NextFacet)[j]=i==j?1.0: 0.0;
        FacetCoords(NextFacet)[DIM]=-coords[i];
        clear_FacetAdj(NextFacet);
        // it is adjacent to the new vertex only
        set_bit(FacetAdj(NextFacet),NextVertex);
        set_bit(VertexAdj(NextVertex),NextFacet);
        for(j=0;j<DIM;j++) if(i!=j){
            set_bit(FacetAdj(NextFacet),j);
            set_bit(VertexAdj(j),NextFacet);
        }
        NextFacet++;
    }
    NextVertex++;
    dd_stats.facetno=DIM; 
}

/* int add_initial_facet(double coords[0:dim])
*     collect all facets of the outer approximation first
*  void add_initial_vertex(int final,double coords[0:dim-1])
*     create a consistent initial polytope from the given vertex.*/
int add_initial_facet(const double coords[/*0 .. DIM */])
{int j; double w;
    if(NextFacet>MaxFacets){
        report(R_fatal,"Resume: more facets than specified\n");
        return 1;
    }
    if(NextVertex>DIM){
        report(R_fatal,"Resume: facets should become before vertices\n");
        return 1;
    }
    w=0.0;
    for(j=0;j<DIM;j++){
        if(coords[j]<0.0){
            report(R_fatal,"Resume: facet %d has negative coefficient\n",NextFacet);
            return 1;
        }
        w += coords[j];
    }
    if(w<PARAMS(PolytopeEps)){
        report(R_fatal,"Resume: facet %d has all zero coefficients\n",NextFacet);
        return 1;
    }
    for(j=0;j<=DIM;j++){
        FacetCoords(NextFacet)[j]=coords[j]/w;
    }
    if(PARAMS(Direction))
        FacetCoords(NextFacet)[DIM] *= -1.0;
    clear_FacetAdj(NextFacet);
    // adjacent to which ideal vertices
    for(j=0;j<DIM;j++){
        if(FacetCoords(NextFacet)[j]<PARAMS(PolytopeEps))
            set_bit(FacetAdj(NextFacet),j);
    }
    NextFacet++;
    dd_stats.facetno++;
    return 0;
}

int add_initial_vertex(int final,const double coords[/* 0 .. DIM-1 */])
{int fno,j; double w;
    if(NextVertex>=MaxVertices){
        report(R_fatal,"Resume: more vertices than specified\n");
        return 1;
    }
    set_in_VertexLiving(NextVertex);
    if(final) set_in_VertexFinal(NextVertex);
    clear_VertexAdj(NextVertex);
    for(fno=0;fno<NextFacet;fno++){ // which facets it is adjacent to
        w=FacetCoords(fno)[DIM];
        for(j=0;j<DIM;j++) w+= coords[j]*FacetCoords(fno)[j];
        if(w<-PARAMS(PolytopeEps)){
            report(R_fatal,"Resume: vertex %d is on the negative side of facet %d (%lg)\n",
               fno,NextVertex-DIM,w);
            return 1;
        }
        if(w<PARAMS(PolytopeEps)){ // adjacent
            set_bit(VertexAdj(NextVertex),fno);
            set_bit(FacetAdj(fno),NextVertex);
        }
    }
    NextVertex++;
    return 0;
}

/* Request memory to accomodate more vertices and facets
*
* void allocate_facet_block()
*    add space for DD_FACET_ADDBLOCK << packshift  more facets; expand
*    FacetCoords, VertexAdj. Adjust FacetBitmapBlockSize
* void allocate_vertex_block(count)
*    add space for 'count' many new vertices
*/

static void allocate_facet_block(void)
{ // extend Facets and VertexBitmap
    dd_stats.facets_allocated_no ++;
    dd_stats.facets_allocated += (DD_FACET_ADDBLOCK<<packshift);
    MaxFacets += (DD_FACET_ADDBLOCK<<packshift);
    FacetBitmapBlockSize += DD_FACET_ADDBLOCK;
    // tell the memmory handling part how much space we need
    yrequest(double,M_FacetCoordStore,MaxFacets,FacetSize);
    yrequest(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize);
    // and do allocation
    if(reallocmem()){ // out of memory
        MaxFacets -= (DD_FACET_ADDBLOCK<<packshift);
        FacetBitmapBlockSize -= DD_FACET_ADDBLOCK;
    }
}

static void allocate_vertex_block(int count)
{int total;
    if(OUT_OF_MEMORY) return;
    for(total=DD_VERTEX_ADDBLOCK<<packshift; total<count;
           total += DD_VERTEX_ADDBLOCK<<packshift);
    dd_stats.vertices_allocated_no ++;
    dd_stats.vertices_allocated += total;
    MaxVertices += total;
    VertexBitmapBlockSize = (MaxVertices+packmask)>>packshift;
    yrequest(double,M_VertexCoordStore,MaxVertices,VertexSize);
    yrequest(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_VertexLiving,1,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_VertexFinal,1,VertexBitmapBlockSize);
    if(reallocmem()){ // out of memory
        VertexBitmapBlockSize = (MaxVertices+packmask)>>packshift;
        MaxVertices -= total;
    }
}

/****************************************************************
* int vertex_intersection(v1,v2)
*    compute the number of facets containing both v1 and v2 */
static int vertex_intersection(int v1,int v2)
{int i,total; register BITMAP_t *L1,*L2;
    total=0; L1=VertexAdj(v1); L2=VertexAdj(v2);
    for(i=0;i<FacetBitmapBlockSize;i++,L1++,L2++){
        total += get_bitcount((*L1)&(*L2));
    }
    return total;
}

/* int get_new_vertexno(int threadID)
*    return next available local vertex number, asking for more memory
*     if necessary. The return value is -1 if cannot allocate memory.
*/
inline static int get_new_vertexno(int threadId) // threadId==0 if ! USETHREADS
{int i;
    i=NewVertex[threadId];
    if(i>=MaxNewVertex[threadId]){ // no more space, ask memory
        MaxNewVertex[threadId] += DD_VERTEX_ADDBLOCK<<packshift;
        trequest(M_NewVertexCoordStore,threadId,MaxNewVertex[threadId]);
        trequest(M_NewVertexAdjStore,threadId,MaxNewVertex[threadId]);
        if(OUT_OF_MEMORY){
           MaxNewVertex[threadId] -= DD_VERTEX_ADDBLOCK<<packshift;
           return -1;
        }
    }
    NewVertex[threadId]++;
    return i;
}

/* double vertex_distance(double *facet, int vno)
*   compute the distance between a hyperplane and a vertex. The
*   hyperplane normal must be non-negative; the first DIM vertices
*   are ideal ones, for those the distance is non-negative.
*/
static double vertex_distance(double *facet,int vno)
{double d=0.0; int i; double *vcoords;
    if(vno<DIM){ // ideal vertex
        return facet[vno];
    }
    vcoords=VertexCoords(vno);
    for(i=0;i<DIM;i++){
        d += (*facet) * (*vcoords);
        facet++; vcoords++;
    }
    d += *facet;
    return d;
}

/* int probe_facet(double *coords)
*     return the score; one with the highest score will be added next.
*     The number of outgoing vertices seems to be a good heuristic */
int probe_facet(double *coords)
{int vno,negvertex;
    dd_stats.probefacet++; negvertex=0;
    for(vno=DIM;vno<NextVertex;vno++) if(is_livingVertex(vno)){
        if(vertex_distance(coords,vno)< - PARAMS(PolytopeEps) ) negvertex++;
    }
    return negvertex;
}

/***********************************************************************
* Recompute vertex coordinates from facets adjacent to it
*
* bool solve_lineq(int d, int DIM+1, double VertexArray[DIM+1,d])
*    solve the system of d linear equation of DIM unknows stored in
*    VertexArray[0:DIM,0:d-1]. Use Gauss elimination: get the largest
*    value in the first column A[*,0], swap this row and the first row;
*    and subtract from all rows so that the first column will be all
*    zero, etc. The matrix rank should be exactly dim (so that we have
*    a non-trivial solution). The solution is returned in A[8..DIM-1].
*    PARAMS(LineqEps) is the rank threshold.
*    Return 1 if the matrix rank is not DIM; otherwise return zero. */

#include "round.h"	/* round the solution to a close rational number */

static int solve_lineq(int d,int DIM1,double *FA)
{int i,jmax,col; double v,vmax;
#define A(i,j)  FA[(i)*DIM1+(j)]
    for(col=0;col<DIM;col++){
       /* get the largest value of column col to A[j,col] */
       jmax=col; vmax=0.0;
       for(i=col;i<d;i++){
           v=A(i,col); if(v<0.0) v=-v;
           if(vmax<v){jmax=i;vmax=v;}
       }
       if(vmax< PARAMS(LineqEps)){ /* too small value */
           report(R_err,"solve_lineq: rank is too small (col=%d, max=%lg)\n"
                 ,col,vmax);
           return 1; /* error */
       }
       if(jmax!=col){ /* swap rows jmax and col */
          for(i=col;i<=DIM;i++){
             v=A(col,i); A(col,i)=A(jmax,i); A(jmax,i)=v;
          }
       }
       v=1.0/A(col,col);
       for(i=col+1;i<=DIM;i++){ A(col,i)*=v; }
       A(col,col)=1.0;
       /* and subtract row col from all rows */
       for(jmax=0;jmax<d;jmax++) if(jmax!=col){
           v=A(jmax,col); A(jmax,col)=0.0;
           for(i=col+1;i<=DIM;i++){ A(jmax,i) -= A(col,i)*v; }
       }
    }
    jmax=-1; vmax=0.0;
    for(i=DIM;i<d;i++){
        v=A(i,DIM);if(v<0.0) v=-v;
        if(v>PARAMS(LineqEps)){ jmax=i; if(vmax<v) vmax=v;}
    }
    if(jmax>0){
        report(R_err,"solve_lineq: rank is too large %lg, increase LineqEps=%lg\n",
                     vmax,PARAMS(LineqEps));
        return 1; /* error */
    }
    /* copy the negated final values to row 0 */
    for(col=0;col<DIM;col++){
        v=-1.0*A(col,DIM); round_to(&v); A(0,col)=v;
    }
    return 0;
#undef A
}

/* void recalucalte_vertex(int info,BITMAP_t *adj,double *old,int threadID)
*    calculate the coordinates of the point which is adjacent to all facets.
*    Complain if out of memory; the system is degenerate, or the old and new
*    coeffs are too far from each other */
static void recalculate_vertex(int info, BITMAP_t *adj,double *old,int threadID)
{double v; int i,j,fno,an,DIM1; BITMAP_t fc; double *FA;
    // the adjacency list has size FacetBitmapBlockSize -- it is facets
    DIM1=DIM+1;
    an=0; for(i=0;i<FacetBitmapBlockSize;i++) an+=get_bitcount(adj[i]);
    // request memory
    talloc2(double,M_FacetArray,threadID,an,DIM1);
    if(OUT_OF_MEMORY) return;
    FA=FacetArray(threadID);an=0;fno=0;
#define A(i,j)	FA[(i)*DIM1+(j)]
    for(i=0;i<FacetBitmapBlockSize;i++){
        j=fno;fc=adj[i]; while(fc){
           while((fc&7)==0){fc>>=3; j+=3;}
           if(fc&1){
               memcpy(&(A(an,0)),FacetCoords(j),DIM1*sizeof(double));
               an++;
           }
           j++; fc>>=1;
        }
        fno += (1<<packshift);
    }
    // all facets are collected, compute their intersection
    if(an<DIM){ // not enough facets, this should never happen
        report(R_err,"recalculate: vertex %d has only %d adjacent facets\n",info,an);
        dd_stats.numerical_error++;
        return;
    }
    if(solve_lineq(an,DIM1,FA)){ // numerical error
        report(R_err,"recalculate: adjacency list of vertex %d is degenerate\n",info);
        dd_stats.numerical_error++;
        return;
    }
    // the new vertex is in A(0,i) for 0<=i<DIM
    for(i=0;i<DIM;i++){
        v=old[i]-A(0,i);
        if(v>PARAMS(VertexRecalcEps) || v<-PARAMS(VertexRecalcEps)){
            report(R_warn,"recalculate: numerical instability at vertex %d,"
               " coord %d (%lg)\n",info,i,v);
            dd_stats.instability_warning++;
        }
        old[i]=A(0,i);
    }
#undef A
}

/* void recalculate_vertices(void)
*    go over all vertices and recalculate their coordinates
*  void thread_recalculate(threadId)
*    split all cases into ThreadNo pieces; each thread executes
*    one of them. If there are no thrads, ThreadNo=1  */

static void thread_recalculate(int threadId) // Id goes from 0 to MaxThreads-1
{int vno,step;
    step=ThreadNo;
    for(vno=DIM+threadId;vno<NextVertex;vno+=step) if(is_livingVertex(vno))
        recalculate_vertex(vno,VertexAdj(vno),VertexCoords(vno),threadId);
}

void recalculate_vertices(void)
#ifdef USETHREADS
{  thread_execute(thread_recalculate); }
#else /* ! USETHREADS */
{   thread_recalculate(0); }
#endif /* USETHREADS */

/**********************************************************************
*
* int is_edge(v1,v2)
*    combinatorial test the check if v1-v2 is and edge. If there are
*    < DIM-1 faces containing both of them, then it is not an edge.
*    Otherwise intersect all facets containing both vertices; it is
*    an edge if the intersection contains no other vertices. The first
*    vertex can be ideal. If no threads, use threadId=0. 
*    This routine assumes that at least 2 facets contain v1 and v2 */
inline static int is_edge(int v1,int v2,int threadId)
{int facetno,i,j,flistlen; BITMAP_t v; BITMAP_t *f0,*f1;
    if(vertex_intersection(v1,v2) < DIM-1)
        return 0; // no  - happens frquently
    /* make a list of all facets adjacent to both v1 and v2,
       and delete v1 and v2 */
    copy_VertexLiving_to(VertexWork(threadId));
    clear_bit(VertexWork(threadId),v1); clear_bit(VertexWork(threadId),v2);
    /* store facets containing both v1 and v2 in VertexList */
    facetno=0; flistlen=0;
    for(i=0;i<FacetBitmapBlockSize;i++){
        if((v=VertexAdj(v1)[i] & VertexAdj(v2)[i])){
            j=facetno;
            while(v){
               while((v&7)==0){ j+=3; v>>=3; }
               if(v&1){
                   FacetList(threadId)[flistlen]=j;
                   flistlen++; 
               }
               j++; v>>=1;
            }
        }
        facetno += (1<<packshift);
    }
    /* now we have all facets adjacent to v1 and v2 in FacetList */
    f0=FacetAdj(FacetList(threadId)[0]); // adjacency list of first facet
    f1=FacetAdj(FacetList(threadId)[1]); // adjacency list of second facet
    for(i=0;i<VertexBitmapBlockSize;i++) if(
       (v=VertexWork(threadId)[i] & f0[i]) && (v &= f1[i])){
         for(j=2;j<flistlen;j++){
           v &= FacetAdj(FacetList(threadId)[j])[i];
         }
         if(v) return 0; // no
    }
    return 1; // yes
}

/* void create_new_vertex(v1,v2)
*    create a new vertex on the edge v1-v2 intersecting the facet ThisFacet
*    Recalculate the vertex coeffs when ExactVertexEq parameter is set */
inline static void create_new_vertex(int v1,int v2,int threadId)
{int newv; double d1,d2,d; int i;
    newv=get_new_vertexno(threadId);
    if(newv<0) return; // no memory
    // adjacency list is the intersection of that of v1 and v2 plus the new facet
    for(i=0;i<FacetBitmapBlockSize;i++)
        NewVertexAdj(threadId,newv)[i] = VertexAdj(v1)[i] & VertexAdj(v2)[i];
    set_bit(NewVertexAdj(threadId,newv),ThisFacet);
    // compute the intersection, v1<0, v2>0
    if(v2<DIM){ // ideal vertex, f[v2]>0
        for(i=0;i<DIM;i++)
            NewVertexCoords(threadId,newv)[i]=VertexCoords(v1)[i];
         d1 = -VertexDist(v1);
         NewVertexCoords(threadId,newv)[v2] += d1/FacetCoords(ThisFacet)[v2];
    } else {    // f(v1)=-d1, f(v2)=d2, (d2*v1+d1*v2)/(d1+d2)
        d1 = -VertexDist(v1); d2 = VertexDist(v2);
        d=1.0/(d1+d2); d1 *= d; d2 *= d;
        for(i=0;i<DIM;i++)
           NewVertexCoords(threadId,newv)[i]=
             d2*VertexCoords(v1)[i]+d1*VertexCoords(v2)[i];
    }
    if(PARAMS(ExactVertex))
        recalculate_vertex(MaxVertices+newv, // report number if error
            NewVertexAdj(threadId,newv),     // adjacency list
            NewVertexCoords(threadId,newv),  // old coordinates, replaced
            threadId);                       // thread
}

/* void make_vertex_living(vno)
*     set the "living" flag for this vertex, add it to the adjacency list
*     of all facets it is adjacent to */
static void make_vertex_living(int vno)
{int fno,i,j; BITMAP_t fc;
    set_in_VertexLiving(vno);
    fno=0;for(i=0;i<FacetBitmapBlockSize;i++){
        j=fno; fc=VertexAdj(vno)[i];
        while(fc){
            while((fc&7)==0){ j+=3; fc>>=3; }
            if(fc&1) set_bit(FacetAdj(j),vno);
            j++; fc>>=1;
        }
        fno += (1<<packshift);
    }
}

/* void search_edges()
*    go over all vertex pairs (v1,v2). v1<0; v2>0 adn check if it is an edge.
*    if yes, add the new vertex, recomputing the coordinates if necessary
*  void thread_search_edges(threadId)
*    split all cases into ThreadNo pieces; each thread executes one of them.
*    If there are no threads, ThreadNo=1, and threadId=0 */

static void thread_search_edges(int threadId) // Id goes from 0 to MaxThreads-1
{int i,j,v1,step; int *PosIdx, *NegIdx;
    step=ThreadNo; // at least 1
    NegIdx=VertexPosnegList+(MaxVertices-1-threadId);
    for(j=threadId;j<dd_stats.vertex_neg;j+=step,NegIdx-=step){
        v1=*NegIdx;PosIdx=VertexPosnegList;
        for(i=0;i<dd_stats.vertex_pos;i++,PosIdx++)
            if(is_edge(v1,*PosIdx,threadId))
                create_new_vertex(v1,*PosIdx,threadId);
    }
}

static void search_edges(void)
#ifdef USETHREADS
{   thread_execute(thread_search_edges); }
#else /* ! USETHREADS */
{   thread_search_edges(0); }
#endif /* USETHREADS */


/* void compress_from(vno)
*    move living vertices from the end to the holes starting from vno */
static void compress_from(int vno)
{int i;
 fill_holes:
    while(vno<NextVertex && is_livingVertex(vno)) vno++;
    while(vno<NextVertex && !is_livingVertex(NextVertex-1)) NextVertex--;
    if(vno<NextVertex){ // vno is empty, NextVertex-1 is used
        NextVertex--;
        move_vertex_to(VertexCoords(NextVertex),VertexAdj(NextVertex),vno);
        make_vertex_living(vno);
        clear_bit(VertexLiving,NextVertex);
        if(is_finalVertex(NextVertex)) set_in_VertexFinal(vno);
        vno++; goto fill_holes;
    }
    // clear VertexFinal
    for(i=0;i<VertexBitmapBlockSize;i++) VertexFinal[i] &= VertexLiving[i];
    // clear adjacency lists again
    for(i=0;i<NextFacet;i++) intersect_FacetAdj_VertexLiving(i);
}

/* void compress_vertices()
*    if the number of vertices decreases, decrease the storage as well */
static void compress_vertices(void)
{int i,vno;
    if(MaxVertices <= get_vertexnum()+2*(DD_VERTEX_ADDBLOCK<<packshift)) return;
    dd_stats.vertex_compressed_no++;
    // clear non-living bits in adjacency list of facets
    for(i=0;i<NextFacet;i++)
        intersect_FacetAdj_VertexLiving(i);
    // find the first free vertex slot
    for(i=vno=0;vno<NextVertex &&(~VertexLiving[i])==0;i++,vno+=(1<<packshift));
    compress_from(vno);
    // compute new limits
    while(MaxVertices > NextVertex+2*(DD_VERTEX_ADDBLOCK<<packshift))
        MaxVertices -= (DD_VERTEX_ADDBLOCK<<packshift);
    VertexBitmapBlockSize = (MaxVertices+packmask)>>packshift;
    yrequest(double,M_VertexCoordStore,MaxVertices,VertexSize);
    yrequest(BITMAP_t,M_VertexAdjStore,MaxVertices,FacetBitmapBlockSize);
    yrequest(BITMAP_t,M_FacetAdjStore,MaxFacets,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_VertexLiving,1,VertexBitmapBlockSize);
    yrequest(BITMAP_t,M_VertexFinal,1,VertexBitmapBlockSize);
    if(reallocmem()){ // fatal error
        report(R_fatal,"Compress vertices error, program aborted\n"); 
        exit(4);
    }
}

/* void request_main_loop_memory(threadID)
*    allocate temporary memory used by a thread in the main loop
*    FacetList:        list of facets adjacent to vertices v1 and v2
*    VertexWork:       bitmap of vertices adjacent to all facets in FacetList
*    NewVertexCoord:   coordinates of new vertices
*    NewVertexAdj:     adjacency bitmap of the a vertex */
inline static void request_main_loop_memory(int threadID)
{   talloc2(int,M_FacetList,threadID,MaxFacets,1);
    talloc2(BITMAP_t,M_VertexWork,threadID,1,VertexBitmapBlockSize);
    talloc2(double,M_NewVertexCoordStore,threadID,
        DD_INITIAL_VERTEXNO,VertexSize);
    talloc2(BITMAP_t,M_NewVertexAdjStore,threadID,
        DD_INITIAL_VERTEXNO,FacetBitmapBlockSize);
    NewVertex[threadID]=0;
    MaxNewVertex[threadID]=DD_INITIAL_VERTEXNO;
}

/* void add_new_facet(double facet[0:DIM])
*    add a new facet which intersects the present approximation. Split
*    vertices into postive, zero, and negative sets relative to their position
*    to the new facet. For each pair of positive/negative vertices check if
*    it is an edge. If yes, add a new vertex at the intersection point. */

/* add a new facet fo the approximation */
void add_new_facet(double *coords)
{double d; int i,j,vno,threadId,AllNewVertex; BITMAP_t fc;
 int *PosIdx, *NegIdx;
    dd_stats.iterations++; dd_stats.facetno++;
    if(NextFacet>=MaxFacets){
        compress_vertices();
        allocate_facet_block();
    }
    // memory for the outer loop
    talloc(double,M_VertexDistStore,MaxVertices,1);
    talloc(int,M_VertexPosnegList,MaxVertices,1);
    for(threadId=0;threadId<ThreadNo;threadId++) request_main_loop_memory(threadId);
    if(OUT_OF_MEMORY){ // indicate that data is till consistent
        dd_stats.data_is_consistent=1;
        return;
    }
    ThisFacet=NextFacet; NextFacet++;
    for(i=0;i<=DIM;i++) FacetCoords(ThisFacet)[i]=coords[i];
    clear_FacetAdj(ThisFacet); // clear the adjacency list
    dd_stats.vertex_pos=0; dd_stats.vertex_neg=0; dd_stats.vertex_zero=0;
    PosIdx = VertexPosnegList; // this goes ahead
    NegIdx = VertexPosnegList+MaxVertices; // this goes backward
    for(vno=0;vno<NextVertex;vno++) if(is_livingVertex(vno)){
       d=VertexDist(vno)=vertex_distance(coords,vno);
       if(d>PARAMS(PolytopeEps)){ // positive size
            *PosIdx=vno; ++PosIdx;
            dd_stats.vertex_pos++;
        } else if(d<-PARAMS(PolytopeEps)){ // negative side
            if(is_finalVertex(vno)){
                report(R_warn,"Final vertex %d is on the negative side of "
                    "facet %d (d=%lg)\n",vno,ThisFacet,d);
                dd_stats.instability_warning++;
// it seems the best thing is to revoke the 'final' flag
                clear_bit(VertexFinal,vno);
            }
            --NegIdx; *NegIdx = vno;
            dd_stats.vertex_neg++;
        } else { // this is adjacent to our new facet
            set_bit(FacetAdj(ThisFacet),vno);
            set_bit(VertexAdj(vno),ThisFacet);
            dd_stats.vertex_zero++;
        }
    }
    if(dd_stats.vertex_neg==0){ // the facet does not cut into the polytope
        dd_stats.vertex_new=0; // no new vertices are added at this step
        if(dd_stats.vertex_zero<DIM){
            report(R_err,"Next facet does not cut into the approximation\n");
            dd_stats.numerical_error++;
            return;
        } else { // check if it is a duplicate
            dd_stats.instability_warning++;
            for(j=ThisFacet-1; j>=0 && j>=ThisFacet-100;j--){
                for(i=0;i<=DIM;i++){
                    d=FacetCoords(j)[i]-coords[i];
                    if(d>PARAMS(PolytopeEps) || d<-PARAMS(PolytopeEps)) break;
                }
                if(i==DIM+1){
                    report(R_warn,"Facet has been added earlier as %d\n",j+1);
                    NextFacet--; // don't show it in the list
                    return;
                }
            } // Not found in the previous 100 additions
            report(R_warn,"Facet %d does not cut into the approximation\n",1+ThisFacet);
        }
        return;
    }
    // some statistics
    d=(double)dd_stats.vertex_pos*(double)dd_stats.vertex_neg;
    if(dd_stats.max_tests<d) dd_stats.max_tests=d;
    dd_stats.avg_tests =((dd_stats.iterations-1)*dd_stats.avg_tests+d)
          /((double)dd_stats.iterations);
    if(DIM <= 2){
        /* search_edges() requires DIIM >=3. For DIM==2, v1-v2 is
           an edge iff there is a facet containing both of them. */
        NegIdx=VertexPosnegList+(MaxVertices-1);
        for(j=0;j<dd_stats.vertex_neg;j++,NegIdx--){
            PosIdx=VertexPosnegList;
            for(i=0;i<dd_stats.vertex_pos;i++,PosIdx++)
                if(vertex_intersection(*NegIdx,*PosIdx)!=0)
                     create_new_vertex(*NegIdx,*PosIdx,0);
        }
    } else { search_edges(); }
    if(OUT_OF_MEMORY || dobreak){
        dd_stats.vertex_new=0;
        return;
    }
    // calculate the number of new vertices to dd_stats.vertex_new
    dd_stats.vertex_new=NewVertex[0];
    for(i=1;i<ThreadNo;i++) dd_stats.vertex_new += NewVertex[i];
    // more statistics
    i=dd_stats.vertex_zero+dd_stats.vertex_pos+dd_stats.vertex_new; // the new vertex number
    if(dd_stats.max_vertices<i) dd_stats.max_vertices=i;
    if(dd_stats.max_vertexadded < dd_stats.vertex_new)
        dd_stats.max_vertexadded=dd_stats.vertex_new;
    dd_stats.avg_vertexadded = ((dd_stats.iterations-1)*dd_stats.avg_vertexadded+
        (double)dd_stats.vertex_new)/((double)dd_stats.iterations);
    // delete negative vertices from VertexLiving
    NegIdx=VertexPosnegList+(MaxVertices-1);
    for(j=0;j<dd_stats.vertex_neg;j++,NegIdx--){ 
         clear_bit(VertexLiving,*NegIdx);
    }
    // move new vertices to NextFacet until there is a space
    AllNewVertex=dd_stats.vertex_new;
    for(threadId=0;threadId<ThreadNo;threadId++){
        while(NewVertex[threadId]>0 && NextVertex<MaxVertices){
            NewVertex[threadId]--; AllNewVertex--;
            move_NewVertex_th(threadId,NextVertex);
            make_vertex_living(NextVertex); NextVertex++;
        }
        if(NextVertex>=MaxVertices) break;
        // threadId remains where we need it later
    }
    if(AllNewVertex==0) return; // done
    while(NewVertex[threadId]==0) threadId++;
    // no more direct space; fill the holes first
    dd_stats.vertex_compressed_no++;
    // clear the adjacency list of facets
    for(i=0;i<NextFacet;i++) intersect_FacetAdj_VertexLiving(i);
    // move new vertices to free vertex slots, if any
    vno=0; for(i=0;i<VertexBitmapBlockSize;i++){
        j=vno; fc=~VertexLiving[i]; // complement ...
        while(fc){
            while((fc&7)==0){ j+=3; fc>>=3; }
            if(fc&1){
                if(j>=MaxVertices) goto finish_compress;
                NewVertex[threadId]--;
                move_NewVertex_th(threadId,j);
                make_vertex_living(j);
                AllNewVertex--;
                if(AllNewVertex==0) goto finish_compress;
                while(NewVertex[threadId]==0) threadId++;
            }
            j++; fc>>=1;
        }
        vno += (1<<packshift);
    }
 finish_compress:
    if(AllNewVertex>0){ // we still have vertices to be included
        allocate_vertex_block(AllNewVertex);
        // if no memory, throw away the rest
        if(OUT_OF_MEMORY) AllNewVertex=0;
        else while(1){
            NewVertex[threadId]--;
            move_NewVertex_th(threadId,NextVertex);
            make_vertex_living(NextVertex);
            NextVertex++;
            AllNewVertex--;
            if(AllNewVertex==0) break;
            while(NewVertex[threadId]==0) threadId++;
        }
        return;
    }
    // until vno all vertex slots are occupied; compress the rest
    compress_from(vno);
}

/***********************************************************************
* int get_next_vertex(from,to[0:DIM-1])
*    return the next living but not final vertex number starting at
*    'from' and store it at 'to'[]. If no such vertex is, return -1. 
*    When from=-1 and RandomVertex is set, pick the starting number
*    randomly. */

/* return a random integer between 0 and v-1 approximately uniformly */
static inline int mrandom(int v)
{ return v<=1 ? 0 : (random()&0x3fffffff)%v; }

int get_next_vertex(int from, double *to/*[0:DIM-1]*/)
{int vno,i,j; BITMAP_t v;
    dd_stats.vertexenquiries++;
    if(from<0){
        if(PARAMS(RandomVertex)){ // start from a random place
            dd_stats.vertexenquiries--;
            vno=get_next_vertex(mrandom(NextVertex),to);
            if(vno>=0) return vno;
        }
        from=0; 
    } else if(from >= NextVertex) return -1;
    j=from & packmask;
    vno = from-j;
    i = vno>>packshift;
    if(j){
        v=(VertexLiving[i] & ~VertexFinal[i])>>j;
        while(v){
            while((v&7)==0){j+=3; v>>=3; }
            if(v&1){
               vno+=j;
               memcpy(to,VertexCoords(vno),DIM*sizeof(double));
               return vno; }
            j++; v>>=1;
        }
        vno += (1<<packshift); i++;
    }
    for(;i<VertexBitmapBlockSize;i++){
        j=0; v= VertexLiving[i] & ~VertexFinal[i];
        while(v){
            while((v&7)==0){j+=3; v>>=3; }
            if(v&1){ 
                vno+=j;
                memcpy(to,VertexCoords(vno),DIM*sizeof(double));
                return vno; }
            j++; v>>=1;
        }
        vno += (1<<packshift);
    }
    return -1; 
}

/***********************************************************************
* Report vertices and facets
*
* void print_vertex(channel,int vno)
*   print vertex 'vno' with prefix 'V' (final) or 'v' (not final)
* void print_vertices(channel)
*   print all vertices
* void print_facet(channel,coords[0:DIM])
*   print scaled and rounded facet coordinates
* void print_facets(channel) 
*   print all facets with prefix 'F' */

#define SCALE_EPS	PARAMS(ScaleEps)

static double round(double x)
/* if x is close to an integer, round it to that integer */
{double res,y;
    if(x<0.0){
        res=(int)(-x+0.5); y=x+res;
        if(y<SCALE_EPS && y>-SCALE_EPS) return -res;
    } else {
        res=(int)(x+0.5); y=x-res;
        if(y<SCALE_EPS && y>-SCALE_EPS) return res;
    }
    return x;
}

static int closetoint(double x)
/* check if x is close to an integer */
{   if(x<0.0) x=-x;
    x -= (int)(x+0.5);
    return (x<SCALE_EPS && x>-SCALE_EPS);
}

static int denum(double x)
/* return d such that d*x is close to an integer */
{int i;
    for(i=1;i<1000;i++)if(closetoint(i*x)) return i;
    return 1;
}

static char *formatvalue(double v)
/* formatted output of a vertex coordinate */
{static char buf[80]; int d;
    round_to(&v); // make it close to a rational, if possible
    if(PARAMS(VertexAsFraction)){
        if(closetoint(v)){sprintf(buf,"%d",(int)(round(v))); return buf; }
        d=denum(v);
        if(d>1){ sprintf(buf,"%d/%d",(int)(round(d*v)),d); return buf; }
    } 
    sprintf(buf,"%.14g",v);
    return buf;
}

void print_vertex(report_type channel,int vno)
/* print vertex vno as a real or as a fraction */
{int j; double dir;
    report(channel,"%c ",is_finalVertex(vno) ? 'V' : 'v');
    dir = PARAMS(Direction) ? -1.0 : 1.0;
    for(j=0;j<DIM;j++){
       report(channel," %s",formatvalue(dir*VertexCoords(vno)[j]));
    }
    report(channel,"\n");
}

void print_vertices(report_type channel)
/* print all vertices */
{int i;
    for(i=DIM;i<NextVertex;i++)if(is_livingVertex(i))
       print_vertex(channel,i);
}

/* int gcd(a,b); int lcm(a,b)
*     compute gcd and lcm of two integers
*  int denom( */
static int gcd(int a,int b) // a,b>=0
{   if(a==0 || a==b) return b;
    if(b==0) return a;
    if(a<b) return gcd(b%a,a);
    return gcd(a%b,b);
}
static int lcm(int a,int b) // a,b>=1
{   if(b==1 || a==b) return a;
    if(a==1) return b;
    return a*(b/gcd(a,b));
}

void print_facet(report_type channel, const double coord[/*0:DIM*/])
/* print a facet using scaling and rounding */
{int j,d; double w,dir;
    dir = PARAMS(Direction) ? -1.0 : 1.0;
    d=1;for(j=0;d<300000 && j<=DIM;j++){ 
        w=coord[j];round_to(&w); d=lcm(d,denum(w));
    }
    for(j=0;j<=DIM;j++){
        w=d*coord[j]; if(j==DIM) w *= dir; round_to(&w);
        report(channel," %.14lg",round(w));
    }
    report(channel,"\n");
}

void print_facets(report_type channel)
{int fno;
    for(fno=0;fno<NextFacet;fno++){
        report(channel,"F ");
        print_facet(channel,FacetCoords(fno));
    }
}

/* void make_checpoint(void)
*     create the next checkpoint file. Fail silently */
void make_checkpoint(void)
{static int version=0;
    open_checkpoint(version);
    report(R_chk,"C checkpoint file #%03d for %s\n",version,PARAMS(ProblemName));
    version++;
    report(R_chk,"C vertices facets rows columns objects\n");
    report(R_chk,"N %d %d %d %d %d\n",NextVertex,get_facetnum()+1,
       PARAMS(ProblemRows),PARAMS(ProblemColumns),PARAMS(ProblemObjects));
    print_facets(R_chk); print_vertices(R_chk);
    close_checkpoint();
}

/* void make_dump(void)
*    when requested by a signal, dump facets and vertices */
void make_dump(void)
{   open_dumpfile();
    report(R_chk,"C snapshot data for %s\n",PARAMS(ProblemName));
    report(R_chk,"C vertices facets rows columns objects\n");
    report(R_chk,"N %d %d %d %d %d\n",NextVertex,get_facetnum()+1,
       PARAMS(ProblemRows),PARAMS(ProblemColumns),PARAMS(ProblemObjects));
    print_facets(R_chk); print_vertices(R_chk);
    close_dumpfile();
}

/************************************************************************
* int check consistency of adjacency bitmaps */

int check_bitmap_consistency(void)
/* check that each vertex has at least DIM adjacent facets
   and each facet has at least DIM adjacent vertices */
{int i,vno,fno,errno; int nn;
    errno=0;
    for(vno=0;vno<NextVertex;vno++) if(is_livingVertex(vno)){
        // number of adjacent facets must be at least DIM
        nn=0;
        for(i=0;i<FacetBitmapBlockSize;i++) nn+=get_bitcount(VertexAdj(vno)[i]);
        if(nn<DIM && (nn<DIM-1 ||vno>=DIM) ){ // except for ideal vertices
           report(R_err,"Vertex %d is on %d facets only (<%d)\n",vno,nn,DIM);
           errno++;
        }
    }
    for(fno=0;fno<NextFacet;fno++){
        nn=0;
        for(i=0;i<VertexBitmapBlockSize;i++) nn+=get_bitcount(FacetAdj(fno)[i]);
        if(nn<DIM){
           report(R_err,"Facet %d contains %d vertices only (<%d)\n",fno,nn,DIM);
           errno++;
        }
    }
    return errno;
}

static int check_vertex_consistency(int vno)
/* check the adjacency bitmaps against the real distance */
{int fno,errno; double d;
    errno=0;
    // go over all vertices
    for(fno=0;fno<NextFacet;fno++){
        d=vertex_distance(FacetCoords(fno),vno);
        if(d<-PARAMS(PolytopeEps)){
            errno++;
            report(R_err,"Vertex %d is on the negative side of facet %d (%lg)\n",
                vno+1,fno+1,d);
        } else if(d<PARAMS(PolytopeEps)){ // adjacent
            if(!extract_bit(FacetAdj(fno),vno)){
                errno++;
                report(R_err,"Facet %d adjacency list: adjacent vertex %d not set\n",
                    fno+1,vno+1);
            }
            if(!extract_bit(VertexAdj(vno),fno)){
                errno++;
                report(R_err,"Vertex %d adjacecny list: adjacent facet %d is not set\n",
                    vno+1,fno+1);
            }
        } else { // not adjacent
            if(extract_bit(FacetAdj(fno),vno)){
                errno++;
                report(R_err,"Facet %d adjacency list: non-adjacent vertex %d is set\n",
                    fno+1,vno+1);
            }
            if(extract_bit(VertexAdj(vno),fno)){
                errno++;
                report(R_err,"Vertex %d adjacency list: non-adjacent facet %d is set\n",
                    vno+1,fno+1);
            }
        }
    }
    return errno;
}

static void thread_check_consistency(int threadId) // 0<=Id<ThreadNo
{int vno,step;
    step=ThreadNo; ErrorNo[threadId]=0;
    for(vno=threadId;vno<NextVertex;vno+=step) if(is_livingVertex(vno))
        ErrorNo[threadId]+=check_vertex_consistency(vno);
}

/** consistency check: adjacency lists are OK **/
int check_consistency(void)
{int i,errno;
  #ifdef USETHREADS
        thread_execute(thread_check_consistency);
  #else /* ! USETHREADS */
        thread_check_consistency(0);
  #endif /* USETHREADS */
    errno=check_bitmap_consistency();
    for(i=0;i<ThreadNo;i++) errno += ErrorNo[i];
    return errno;
}


/* EOF */


