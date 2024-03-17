/** outer.c outer approximation using glpk oracle **/

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include "report.h"
#include "outer.h"
#include "data.h"
#include "params.h"
#include "poly.h"
#include "glp_oracle.h"

/*********************************************************************
* Miscellaneous: initialize random; get elapsed time; nice printing
*
* void initialize_random()
*    initialize random by asking the time and the pid of the process.
*    The random() function returns an integer in the range [0:2^{31}-1]
*
* unsigned long timenow
*    elapsed time since the first call of gettime100()
*
* unsigned long gettime100()
*    set timenow and return the elapsed time
*
* char *showtime(unsigned long t)
*    readable printout of the time in days, hours, minutes and seconds
*
* char *readable(double w, int slot)
*    convert w to a string with k,M,G,T postfix using one of the 
*    static strings at slots 0,1,2,3
*/

#include <sys/time.h> /* gettimeofday() */
#include <unistd.h>   /* pid() */
/** initialize random **/
inline static void initialize_random(void)
{struct timeval tv;
    gettimeofday(&tv,NULL);
    srandom( PARAMS(TrueRandom) ? getpid()^tv.tv_sec^tv.tv_usec : 0x12345678 );
}
/** the elapsed time in 0.01 seconds */
unsigned long timenow=0;
/** set timenow and return the elapsed time */
static unsigned long gettime100(void)
{struct timeval tv; static unsigned long starttime=0;
    if(gettimeofday(&tv,NULL)) return timenow; // some problem
    if(starttime){
        timenow = (tv.tv_sec*100 + (tv.tv_usec+5000u)/10000u)-starttime;
        return timenow;
    }
    starttime=tv.tv_sec*100 + (tv.tv_usec+5000u)/10000u;
    timenow = 0;
    return timenow;
}
/** print time in pretty format **/
static char *showtime(unsigned long t)
{static char buff[50]; int m,s;
    if(t<60*100){ sprintf(buff,"%.2f",((double)t)*0.01); return buff;}
    s=(int)((t+50)/100); m=s/60; s=s%60;
    if(m<60){ sprintf(buff,"%d:%02d",m,s); return buff; }
    if(m<60*24){ sprintf(buff,"%d:%02d:%02d",m/60,m%60,s); return buff;}
    sprintf(buff,"%dd%02d:%02d:%02d",m/1440,(m/60)%24,m%60,s); return buff;
}
/** print a huge double with postfix k,M,G,T to a slot 0..3 **/
static char* readable(double w, int slot)
{static char slots[4][30]; char *buff;
    buff=slots[slot];
    if(w<0.0){ w=0.0; } /* only >=0 numbers */
    if(w<1000.0){ sprintf(buff,"%.2lf",w); return buff;}
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfk",w); return buff; }
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfM",w); return buff; }
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfG",w); return buff; }
    w*=0.001;
    if(w<1000.0){ sprintf(buff,"%.2lfT",w); return buff; }
    /* too large value */
    sprintf(buff,"%lgT",w); return buff;
}

/*******************************************************************
* Progress report
*
* unsigned long progresstime, progressdelay
*    when was the last progress report made
*    minimum delay between two reports, in 0.01 seconds
*
* unsigned long chktime, chkdelay
*    when was the last checkpoint created, minimum delay between
*    two checkpoints
*
* int poolstat
*    number of entries in FacetPool, set by find_next_facet()
*
* int vertexstat
*    whether vertex statistics changed calling add_new_vertex()
*
* void progress_stat(void)
*    show progress report, save the report time
*
* void progress_stat_if_expired(int ask_time)
*    show progress report when expired
*
* void report_new_facet(int ask_time)
*    report coordinates of the new facet when requested; give 
*    status report when the delay has passed
*
* void report_new_vertex(int vertexno)
*    report vertices for the new final vertex; give status 
*    report when the delay has passed
*
* void report_memory(void)
*    when total memory allocated changes, call report_memory_usage()
*
* int limit_reached(void)
*    return 1 if the allocated memory exceeds the one set in PARAMS(MemoryLimit)
*    or the time limit exceeds the one set in PARAM(TimeLimit)
*/
static unsigned long progresstime=0, progressdelay=0;
static unsigned long chktime=0, chkdelay=0;

static int poolstat=0, vertexstat=0;

static void progress_stat(void)
{   progresstime=timenow;
    get_dd_vertexno(); // fills living_vertex_no and final_vertex_no
    report(R_txt,"I%8.2f] Elapsed: %s, facets: %d, vertex final: %d, pending: %d",
        0.01*(double)timenow,
        showtime(timenow),
        dd_stats.facetno,
        dd_stats.final_vertex_no,
        dd_stats.living_vertex_no-dd_stats.final_vertex_no);
    if(poolstat){ report(R_txt,", pool: %d",poolstat); poolstat=0; }
    if(vertexstat){
        report(R_info,", eq: %d, out: %d, in: %d",
           dd_stats.vertex_zero, dd_stats.vertex_neg,dd_stats.vertex_new);
        vertexstat=0;
    }
    report(R_txt,"\n");
}

static void progress_stat_if_expired(int ask_time)
{   if(PARAMS(ProgressReport)==0) return;
    if(ask_time) gettime100();
    if(timenow-progresstime < progressdelay) return;
    progress_stat();
    flush_report();
}

static void report_new_facet(int ask_time)
{int flush=0;
    if(PARAMS(ProgressReport)==0 && !PARAMS(FacetReport)) return;
    if(PARAMS(ProgressReport)){
      if(ask_time){ gettime100(); ask_time=0; }
      if(timenow-progresstime >= progressdelay){
        progress_stat();
        flush=1;
      }
    }
    if(PARAMS(FacetReport)){
      if(ask_time) gettime100();
      report(R_txt,"[%8.2f] F ",0.01*(double)timenow);
      print_facet(R_txt,OracleData.ofacet);
      flush=1;
    }
    if(flush) flush_report();
}

static void report_new_vertex(int vno)
{int flush=0;
    if(PARAMS(ProgressReport)==0 && PARAMS(VertexReport)==0) return;
    gettime100();
    if(PARAMS(ProgressReport)){
       if(timenow-progresstime >= progressdelay){
         progress_stat(); flush=1;
       }
    }
    if(PARAMS(VertexReport)){
      report(R_txt,"[%8.2f] ",0.01*(double)timenow);
      print_vertex(R_txt,vno); flush=1;
    }
    if(flush) flush_report();
}

static void report_memory(void)
{static int last_memreport=0; char buff[50];
    if( PARAMS(MemoryReport)<2 ||
       last_memreport==dd_stats.memory_allocated_no ||
       dd_stats.out_of_memory ) return;
    last_memreport=dd_stats.memory_allocated_no;
    sprintf(buff,"I%8.2f] Memory allocation table", 0.01*(double)timenow);
    report_memory_usage(R_txt,0,buff);
    flush_report();
}

inline static int limit_reached(void)
{  if(PARAMS(TimeLimit)>=60 && timenow > 100ul*(unsigned long)PARAMS(TimeLimit))
      return 1;
   if(PARAMS(MemoryLimit)>=100 && 
      (((double)dd_stats.total_memory)*1e-6 > ((double)PARAMS(MemoryLimit))))
      return 1;
   return 0;
}

/*********************************************************************
* Print statistics
*
* EQSEP, DASHSEP
*   separators made of = and -
*
* void dump_and_save(int status)
*   dump and save vertices and facets. The termination status:
*     0  - normal
*     1  - error but data is consistent
*     2  - error with inconsistent data
*     3  - interrupt
*/

#define EQSEP   "================================"
#define DASHSEP "--------------------------------"

static void dump_and_save(int status)
{unsigned long endtime; int partial;
    endtime=gettime100(); // program finished
    if(PARAMS(ProgressReport)) progress_stat();
    partial= status==0 ? 0 : 1; // print data when completed
    if(PARAMS(PrintVertices) > partial){
        if(partial) report(R_txt,"Partial list of vertices:\n");
        report(R_txt,"\n" DASHSEP "\n");
        print_vertices(R_txt);
    }
    if(PARAMS(PrintFacets) > partial){
        if(partial) report(R_txt,"Partial list of facets:\n");
        report(R_txt,"\n" DASHSEP "\n");
        print_facets(R_txt);
    }
    if(PARAMS(PrintStatistics)){ // statistics, only if not quiet
      int oraclecalls,oraclerounds; unsigned long oracletime;
      get_oracle_stat(&oraclecalls,&oraclerounds,&oracletime);
      report(R_txt,"\n" DASHSEP "\n"
      "Problem %s\n"
      " name                    %s\n"
      " output                  %s\n"
      "%s%s%s"
      "%s%s%s"
      " rows, cols, objs        %d, %d, %d\n"
      " vertices, facets        %d, %d\n",
      status==0 ? "completed" : status<=2 ? "aborted with error" : "interrupted",
      PARAMS(ProblemName),
      PARAMS(SaveFile) ? PARAMS(SaveFile) : 
        PARAMS(SaveVertexFile)||PARAMS(SaveFacetFile) ? "" : "[none]",
      PARAMS(SaveVertexFile) ? "   vertices              " : "",
      PARAMS(SaveVertexFile) ? PARAMS(SaveVertexFile) : "",
      PARAMS(SaveVertexFile) ? "\n": "",
      PARAMS(SaveFacetFile) ?  "   facets                " : "", 
      PARAMS(SaveFacetFile) ? PARAMS(SaveFacetFile) : "",
      PARAMS(SaveFacetFile) ? "\n" : "",
      PARAMS(ProblemRows), PARAMS(ProblemColumns), PARAMS(ProblemObjects),
      get_vertexnum(), get_facetnum());
      report(R_txt, " total time              %s\n",
         showtime(endtime)); 
      report(R_txt, DASHSEP "\nStatistics\n"
      "LP\n"
      " oracle calls            %d\n"
      "   avg iterations/call   %s\n"
      " total oracle time       %s\n",
      oraclecalls, readable((0.0001+oraclerounds)/(0.0001+oraclecalls),0),
      showtime(oracletime));
      report(R_txt,
      "Combinatorics\n"
      " vertices probed         %d\n"
      " vertices # max          %d\n"
      " facets added            %d\n"
      " memory allocated        %s%s\n",
      dd_stats.vertexenquiries, dd_stats.max_vertices,
      dd_stats.iterations+1,
      readable(dd_stats.max_memory,0),
      dd_stats.out_of_memory ? " (out of memory)" : "");
#ifdef USETHREADS
      report(R_txt, " threads                 %d\n",PARAMS(Threads));
#endif      
      if(dd_stats.instability_warning) report(R_txt,
      " instability warnings    %d\n",
      dd_stats.instability_warning);
      report(R_txt,
      " storage expansion\n"
      "   vertices # / upto     %d / %d\n"
      "   facets   # / upto     %d / %d\n"
      " vertex manipulating\n"
      "   added avg / max       %s / %s\n"
      "   compressed            %d\n"
      " number of edge tests\n"
      "   avg / max             %s / %s\n",
      dd_stats.vertices_allocated_no,dd_stats.vertices_allocated,
      dd_stats.facets_allocated_no,dd_stats.facets_allocated,
      readable(dd_stats.avg_vertexadded,0),readable(dd_stats.max_vertexadded,1),
      dd_stats.vertex_compressed_no,
      readable(dd_stats.avg_tests,2),readable(dd_stats.max_tests,3));
      if(PARAMS(MemoryReport)>0 || dd_stats.out_of_memory)
         report_memory_usage(R_txt,1,"Memory allocation:");
      if(PARAMS(PrintParams))
         show_parameters("Parameters with non-default values:\n");
    } else {
        if(PARAMS(MemoryReport)>0 || dd_stats.out_of_memory)
            report_memory_usage(R_txt,1,"\n" DASHSEP "\nMemory allocation:");
        if(PARAMS(PrintParams))
            show_parameters(DASHSEP "\nParameters with non-default values:\n");
    }
    partial = status<2 ? 0 : 1; // save data when consistent
    if(PARAMS(SaveFacets)>partial || (partial==0 && PARAMS(SaveFacetFile))){
        report(R_savefacet,"C name=%s, rows=%d, cols=%d, objs=%d\n"
          "C vertices=%d, facets=%d\n\n",
          PARAMS(ProblemName), 
          PARAMS(ProblemRows), PARAMS(ProblemColumns), PARAMS(ProblemObjects),
          get_vertexnum(), get_facetnum());
        if(status) report(R_savefacet,"C *** Partial list of facets ***\n");
        print_facets(R_savefacet);
        report(R_savefacet,"\n");
    }
    if(PARAMS(SaveVertices)>partial || (partial==0 && PARAMS(SaveVertexFile))){
        report(R_savevertex,"C name=%s, rows=%d, cols=%d, objs=%d\n"
          "C vertices=%d, facets=%d\n\n",
          PARAMS(ProblemName), 
          PARAMS(ProblemRows), PARAMS(ProblemColumns), PARAMS(ProblemObjects),
          get_vertexnum(), get_facetnum());
        if(status) report(R_savevertex,"C *** Partial list of vertices ***\n");
        print_vertices(R_savevertex);
        report(R_savevertex,"\n");
    }
    close_savefiles();
}

/*********************************************************************
* Find the next facet to be added
*
* int FacetPoolAfter = 20
* int FacetPoolMinVertices = 1000
*   use facetpool only when we have at least that many facets, or
*   that many unprocessed vertices
*   
* facetpool_t facetpool[FacetPoolSize]
*   facets known but not added yet to the approximation.
*
* int vertices_recalculated
*   TRUE just after calling recalculate_vertices(); and set to FALSE
*   when calling add_new_vertex()
*
* int init_facetpool()
*   allocates memory to the facet pool. Return non-zero if out of
*   memory.
*
* int same_vector(int dim, double v1[0:dim-1], double v2[0:dim-1])
*   checks whether the two vertices are the same (1) or not (0)
*
* int next_facet_coords(int checkFacetPool)
*   try to generate a new facet. Ask the oracle about the next
*   vertex returned by get_next_vertex(-1). Check whether the vertex has
*   been asked before if checkFacetPool is set. If the returned
*   facet is on the vertex, mark the vertex as final and repeat.
*   Return value:
*     0:  there are no more vertices, the algorithm finished
*     1:  interrupted
*     4:  some error (oracle failed, computational error, etc)
*     5:  next (random) vertex is already in FacetPool (only
*            if checkFacetPool!=0)
*     6:  next facet is in OracleData.ofacet (maybe from bootfile)
*     7:  memory or time limit exceeded
*
* int fill_facetpool(int limit)
*   fill the facet pool, limiting unsuccessful oracle calls to limit.
*   Return value:
*     0:  pool is filled (maybe empty if no more vertex)
*     1:  interrupt
*     4:  numerical error
*     7:  memory or time limit exceeded
*
* int find_next_facet(void)
*   find the next facet to be added to the approximating polytope.
*   Without facet pool it calls next_facet_coords().
*   Otherwise returns the facet with the largest probe_facet(v)
*   score. Return values are the same as for next_facet_coords().
*/

#define FacetPoolAfter 20  /* use vertexpool only after that many facets */
#define FacetPoolMinVertices 1000 /* or after that many unprocessed vertices */

typedef struct {
    int    occupied;    /* 0=no, 1=yes */
    double *vertex;     /* pointer to vertex which was asked */
    double *facet;      /* pointer to the facet eq */
} facetpool_t;

static facetpool_t *facetpool=NULL;

#define DIM     PARAMS(ProblemObjects) /* problem dimension */

static int vertices_recalculated=0;  /* set after vertices are recalculated */

static int init_facetpool(void) /* call only when DIM has been set */
{int i; double *pool;
    if(PARAMS(FacetPoolSize)<5) return 0; // don't use it
    facetpool=malloc(PARAMS(FacetPoolSize)*sizeof(facetpool_t));
    pool=malloc(PARAMS(FacetPoolSize)*(2*DIM+1)*sizeof(double));
    if(!facetpool || !pool){
        report(R_fatal,"init_facetpool: out of memory\n");
        return 1;
    }
    for(i=0;i<PARAMS(FacetPoolSize);i++){
        facetpool[i].occupied=0;
        facetpool[i].vertex=pool; pool+=DIM;   // question
        facetpool[i].facet=pool;  pool+=DIM+1; // answer
    }
    return 0;
}

inline static int same_vector(int dim, const double f1[], const double f2[])
{int i; double d;
    for(i=0;i<dim;i++){
       d=f1[i]-f2[i];
        if(d> PARAMS(PolytopeEps) || d < -PARAMS(PolytopeEps) )
            return 0; /* no */
    }
    return 1; /* yes */
}

static int next_facet_coords(int checkFacetPool)
{int i,j; double d; int boottype; int vertex_err=0;
again:
    if(dobreak) return 1; /* interrupt meanwhile */
    if(limit_reached()) return 7; /* memory or time limit */
    // boot file
    while(nextline(&boottype)) if(boottype==3){ // F line
        if(parseline(DIM+1,OracleData.ofacet)){
            return 4; /* error */
        }
        memset(OracleData.overtex,0,DIM*sizeof(double));
/* this makes an oracle reply to the all zero vertex query when filling
   the facetpool. Thus if get_next_vertex() returns the all zero point
   and checkFacetPool=1, then returns 5, and then this query will not
   be asked from the oracle. This happens until all such entries are
   purged from the pool. */
        return 6; /* facet is OK */
    }
    j=get_next_vertex(-1,OracleData.overtex);
    if(j<0) return 0; /* terminated successfully */
    if(checkFacetPool){ /* ask oracle only when not asked before */
        for(i=0;i<PARAMS(FacetPoolSize);i++) if(facetpool[i].occupied
           && same_vector(DIM,facetpool[i].vertex,OracleData.overtex)){
            return 5; // vertex in OracleData.overtex was asaked before
        }
    }
    i=ask_oracle();
    if(i==ORACLE_UNBND){ // on the boundary
        mark_vertex_as_final(j);
        report_new_vertex(j);// progress_stat_if_expired(1);
        goto again;
    }
    if(i!=ORACLE_OK){ // oracle returned with an error
        return  4;    // oracle failed
    }
    // OracleData.ofacet is normalized facet eq
    d=OracleData.ofacet[DIM];
    for(i=0;i<DIM;i++){
        d+= OracleData.ofacet[i]*OracleData.overtex[i];
    }
    if(d > PARAMS(PolytopeEps)){ /* numerical error */
        // make sure vertices are recalculated immediately before issuing an error
        if(vertices_recalculated){
            vertex_err++;
            if(vertex_err<3){
                report(R_warn,"Vertex %d is inside the polytope (d=%lg), trying "
                     "another vertex ...\n",j,d);
                dd_stats.instability_warning++;
                goto again;
            }
            report(R_fatal,"Vertex %d is inside the polytope (d=%lg)\n",j,d);
            return 4;
        }
        report(R_warn,"Vertex %d is inside the polytope, d=%lg, recalculating "
                      "vertices ...\n",j,d);
        dd_stats.instability_warning++;
        recalculate_vertices();
        if(dd_stats.out_of_memory || dd_stats.numerical_error){
            return 4; // error during computation
        }
        vertices_recalculated=1;
        goto again;
    }
    if(d>-PARAMS(PolytopeEps)){ /* vertex is on the facet */
        if(d<0) d=-d;
        report(R_warn,"Numerical instability at vertex %d, distance=%lg < PolytopeEps\n",
           j,d);
        dd_stats.instability_warning++;
        mark_vertex_as_final(j);
        report_new_vertex(j);// progress_stat_if_expired(1);
        goto again;
    }
    return 6;
}

static int fill_facetpool(int limit)
{int i,ii; int oracle_calls=0;
    for(i=0;i<PARAMS(FacetPoolSize);i++)if(!facetpool[i].occupied)
        switch(next_facet_coords(1)){
      case 0:  return 0; /* no more vertices or done */
      case 1:  return 1; /* break */
      case 4:  return 4; /* error */
      case 7:  return 7; /* memory or time limit */
      case 5:  break;    /* the vertex has been encoutered again, skip */
      default:           /* 6, facet is in OracleData */
            for(ii=0;ii<PARAMS(FacetPoolSize);ii++) if(facetpool[ii].occupied
               && same_vector(DIM,OracleData.ofacet,facetpool[ii].facet)) break;
            if(ii==PARAMS(FacetPoolSize)){ // the facet is not in FacetPool
                memcpy(facetpool[i].vertex,OracleData.overtex,DIM*sizeof(double));
                memcpy(facetpool[i].facet,OracleData.ofacet,(DIM+1)*sizeof(double));
                facetpool[i].occupied=1;
            } else { // got the same facet, change vertex to the latter one
                memcpy(facetpool[ii].vertex,OracleData.overtex,DIM*sizeof(double));
                memcpy(facetpool[ii].facet,OracleData.ofacet,(DIM+1)*sizeof(double));
                // and check how many unsuccessful calls were made
                oracle_calls++;
                if(limit && oracle_calls>=limit) return 0; // done
            }
    }
    return 0;
}

static int find_next_facet(void)
{int i,maxi,cnt; int w,maxw;
    if(PARAMS(FacetPoolSize)<5 || ( dd_stats.facetno<FacetPoolAfter &&
      dd_stats.vertex_zero+dd_stats.vertex_pos+dd_stats.vertex_new<FacetPoolMinVertices))
       return next_facet_coords(0);
    // fill the facet pool
    i=fill_facetpool(PARAMS(OracleCallLimit));
    if(i) return i; // some error
    /* find the score of stored facets */
    maxi=-1; maxw=0; cnt=0;
    for(i=0;i<PARAMS(FacetPoolSize);i++) if(facetpool[i].occupied){
        cnt++;
        w=probe_facet(facetpool[i].facet);
        if(maxi<0 || maxw<w){ maxi=i; maxw=w; }
    }
    poolstat=cnt;
    if(maxi<0) return 0; // no more vertices
    facetpool[maxi].occupied=0;
    memcpy(OracleData.ofacet,facetpool[maxi].facet,(DIM+1)*sizeof(double));
    return 6; // next facet is in OracleData.ofacet
} 

/**************************************************************************
* The main loop of the algorithm
* int outer(void)
*   when it starts, all parameters in PARAMS have been set. The steps are
*   o  load_vlp() reads in the the MOLP problem form a vlp file
*   o  check that output files are writable
*   o  initilize_oracle_() sets the oracle parameters, check feasibility
*   o  the first approximation comes from two sources. The "resume"
*        file contains facets (first) and non-ideal vertices (those marked
*        'final' are not queried from the oracle) of the initial polytope.
*        call init_dd_structure(vertices,facets), then add_initial_facet(),
*        add_initial_vertex(), finally check_bitmap_consistency()
*        In case of no 'resume' file the oracle is asked to provide an
*        external point; call init_dd(vertex) to set up the polytope with
*        the positive ideal endpoints of the coordinate axes.
*   loop:
*   o  get_next_vertex(-1) returns the next vertex of the approximating
*        polytope which is not known to be final. If none, the algorithm
*        terminates. Otherwise
*   o  ask_oracle() returns a facet separating the vertex from the final
*        polytope. If the vertex is on the facet, then
*   o  mark the vertex as final, and goto loop.
*        If the vertex is on the negative side of the facet then
*   o  add_new_facet(), and goto loop.
*   Return values:
*     0:  done
*     1:  data error before start (no input/output file, syntax error)
*     2:  problem is unbounded (only at the first approximation)
*     3:  no feasible solution
*     4:  computational error (filling vertex pool, initial data, out of memory)
*     5-7: from break_outer()
*
* void report_error(void)
*   report error fatal error
*
* int break_outer(int how)
*   when the outer routine is interrupted by the signal, this procedure
*   kicks in. It goes over all vertices of the actual approximation,
*   and calls the oracle whether it is final or not. The argument
*   tells if interrupt or memory limit reached.
*     5:  normal termination: no postprocessing necessary or done
*     6:  error during postprocess
*     7:  postprocessing aborted
*
* int handle_new_facet(void)
*   the new facet is in OracleData.ofacet; make reports, add as a new 
*   facet by calling add_new_facet(), take care of timed actions such as
*   reports and checkpoints.
*   Return values:
*     1:  OK
*     0:  error during computation (cannot continue)
*/

/** error report **/
static void report_error(void)
{
    gettime100();
    report(R_fatal,"\n\n" EQSEP "\n%s after %s, vertices: %d, facets: %d\n",
       dd_stats.out_of_memory ? "Out of memory" : "Numerical error",
       showtime(timenow),get_vertexnum(), get_facetnum());
    if(dd_stats.out_of_memory)
       report_memory_usage(R_fatal,1,"Memory allocation table");
}

/* the main loop was interrupted; returns 5,6,7 */
static int break_outer(int status)
{int i,j; unsigned long aborttime;
    aborttime=gettime100(); dobreak=0;
    if(status>0 && PARAMS(TimeLimit)>=60 &&
        aborttime > 100ul*(unsigned long)PARAMS(TimeLimit)) status=2;
    report(R_fatal,"\n\n" EQSEP "\n%s after %s, vertices: %d, facets: %d\n",
       status==0 ? "Program run was interrupted" :
       status==1 ? "Memory limit reached" : "Time limit reached",
       showtime(aborttime), get_vertexnum(), get_facetnum());
    if(!PARAMS(ExtractAfterBreak)) return 5; // normal termination
    if(PARAMS(PrintVertices)<2 && PARAMS(SaveVertices)<2 &&
       !PARAMS(VertexReport)){
        report(R_fatal,"Result of postpricessing would be lost, not doing...\n");
        return 5;
    }
    report(R_fatal,"Checking additional vertices. This may take some time...\n"
        EQSEP "\n\n");
    // at least one of PrintVertices and SaveVertices should be set
    if(PARAMS(PrintVertices)<2 && PARAMS(SaveVertices)<2)
        PARAMS(PrintVertices)=2;
    free_adjacency_lists(); /* release adjacency lists */
    j=-1;
    while((j=get_next_vertex(j+1,OracleData.overtex))>=0){
        if(dobreak){
            dobreak=0;
            if(gettime100()-aborttime>50){
                report(R_fatal,"\n" EQSEP "\n"
                  "Post-processing aborted after %s\n",
                  showtime(timenow-aborttime));
                return 7; // postprocess aborted
            }
        }
        i=ask_oracle();
        if(i==ORACLE_UNBND){ // on the boundary
            mark_vertex_as_final(j);
            report_new_vertex(j);
        } else if(i!=ORACLE_OK) {
            return 6; // error during postprocess
        } else { // vertex is out
            progress_stat_if_expired(1);
        }
    }
    return 5; // terminated
}

static int handle_new_facet(void)
{   report_new_facet(0);  // progress report
    add_new_facet(OracleData.ofacet);
    gettime100(); vertexstat=1; progress_stat_if_expired(0);
    if(dd_stats.out_of_memory || dd_stats.numerical_error)
        return 0; // error meanwhile
    vertices_recalculated=0;
    // recalculate if instructed so
    if(PARAMS(RecalculateVertices)>=5 &&
       ((1+dd_stats.iterations)%PARAMS(RecalculateVertices))==0){
        report(R_info,"I%8.2f] recalculating vertices...\n",0.01*(double)timenow);
        recalculate_vertices();
        gettime100();
        if(dd_stats.out_of_memory || dd_stats.numerical_error){
            return 0; // error during computation
        }
        vertices_recalculated=1;
    }
    if(PARAMS(CheckConsistency)>=5 &&
       ((1+dd_stats.iterations)%PARAMS(CheckConsistency))==0){
        // report what we are going to do
        report(R_warn,"I%8.2f] checking data consistency...\n",
              0.01*(double)timenow);
        if(check_consistency()) { // error
            report(R_fatal,"Consistency error: data structure has numerical errors.\n");
            return 0;
        }
        gettime100();
    }
    report_memory();
    return 1; // OK
}

typedef enum {	/* where data come from */
  inp_none,	/* no special input */
  inp_boot,	/* reading facets from --boot */
  inp_resume	/* facets and vertices from --resume */
} input_type_t;

int outer(void)
{int retvalue=0; input_type_t inp_type=inp_none;
    initialize_random();  // initialize random numbers
    if(load_vlp()) return 1; // data error before start
    if(check_outfiles()) return 1;
    if(PARAMS(BootFile)){ // we have a bootfile
        if(init_reading(PARAMS(BootFile))) return 1; 
        inp_type=inp_boot;
    } else if(PARAMS(ResumeFile)){
        if(init_reading(PARAMS(ResumeFile))) return 1;
        inp_type=inp_resume;
    }
    report(R_info,"C MOLP problem=%s, %s\n"
        "C rows=%d, columns=%d, objectives=%d\n",
        PARAMS(ProblemName),PARAMS(Direction)?"maximize":"minimize",
        PARAMS(ProblemRows),PARAMS(ProblemColumns),PARAMS(ProblemObjects));
    gettime100();  // initialize elapsed time
    if(init_facetpool()) return 1;
    progressdelay = 100*(unsigned long)PARAMS(ProgressReport);
    switch(initialize_oracle()){  // initialize oracle
      case ORACLE_OK:	break;	  // OK
      case ORACLE_EMPTY:return 3; // no feasible solution
      default:		return 4; // oracle error, message given
    }
    chkdelay = 100*(unsigned long)PARAMS(CheckPoint);
    if(inp_type==inp_resume){ // read initial polytope from ResumeFile
        int linetype=0; double args[5];
        if(!nextline(&linetype) || linetype!=4 || parseline(5,&args[0])){
            report(R_fatal,"Resume: missing N line in file '%s'\n",
                PARAMS(ResumeFile));
            return 1; 
        }
        // rows, columnd objects
        if(args[2]!=(double)PARAMS(ProblemRows) ||
           args[3]!=(double)PARAMS(ProblemColumns) ||
           args[4]!=(double)PARAMS(ProblemObjects)) {
            report(R_fatal,"Resume: file '%s' belongs to a different MOLP problem\n",
                 PARAMS(ResumeFile));
            return 1; 
        }
        if(init_dd_structure((int)args[0],(int)args[1])) return 1;
        // and read initial approximation
        while(nextline(&linetype)) switch(linetype){
           case 1:	// 'V' line, after 'F' lines
             if(parseline(DIM,OracleData.overtex) ||
                 add_initial_vertex(1,OracleData.overtex)) return 1;
             break;
           case 2:      // 'v' line, after 'F' lines
             if(parseline(DIM,OracleData.overtex) ||
                 add_initial_vertex(0,OracleData.overtex)) return 1;
             break;
           case 3:	// 'F' line, these should come first
             if(parseline(DIM+1,OracleData.ofacet) ||
                 add_initial_facet(OracleData.ofacet)) return 1;
             break;
           default:	// ignore
             break;
        }
        if(check_bitmap_consistency()){ // in particular complains if no data
            report(R_fatal,"Resume: incosistent data in file %s\n",
                PARAMS(ResumeFile));
            return 1;
        }
    } else {  // create the initial approximation
        if(init_dd_structure(0,0)) return 1;
        switch(get_initial_vertex()){
           case ORACLE_OK:   break;    // OK
           case ORACLE_UNBND:return 2; // problem unbounded, message given
           default:          return 4; // oracle error
        }
        init_dd(OracleData.overtex);  // create the first approximation
    }
#ifdef USETHREADS
    if(create_threads()) return 1;
#endif
    chktime=gettime100(); // last checktime
    progress_stat_if_expired(0);
again:
    gettime100();
    if(dodump){ // request for dump
        dodump=0;
        report(R_info,"I%08.2f] Dumping vertices and facets\n",0.01*(double)timenow);
        make_dump();
    }
    if(PARAMS(CheckPointStub) && chktime+chkdelay<=timenow){
        chktime=timenow;
        report(R_info,"I%08.2f] Checkpoint: dumping vertices and facets\n",
            0.01*(double)chktime);
        make_checkpoint();
    }
    switch(find_next_facet()){
      case 0: /* no more facets */
        dump_and_save(0); retvalue=0; goto leave;
      case 1: /* interrupt */
        retvalue=break_outer(0);
        dump_and_save(3); goto leave;
      case 4: /* numerical or other error */
        report_error();
        dump_and_save(2); retvalue=4; goto leave;
      case 7: /* memory or time limit */
        retvalue=break_outer(1);
        dump_and_save(3); goto leave;
      default: /* next facet returned */
        if(handle_new_facet()) goto again;
        dump_and_save(dd_stats.data_is_consistent? 1 : 2);
        retvalue=4; goto leave;
    }
leave:
#ifdef USETHREADS
    stop_threads();
#endif
    return retvalue;
}

/* EOF */
