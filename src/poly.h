/** poly.h --  combinatorial part using double description method **/

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

/**********************************************************************
* Constant values
*
* BITMAP_64, BITMAP_32
*    bitmaps use either 64 (default) or 32 bit wide words. Use
*    -DBITMAP32 as compiler argument to use 32 bits in bitmaps.
*
* int MAXIMAL_ALLOWED_DIMENSION
*    the maximal dimension we are willing to handle.
*
* int DD_INITIAL_VERTEXNO
* int DD_INITIAL_FACETNO
*    initial number of vertices and facets.
*
* int DD_VERTEX_ADDBLOCK
* int DD_FACET_ADDBLOCK
*    space for vertices and facets are allocated in large blocks.
*    These values multiplied by bitmap word size (64 or 32) determine
*    how many new vertices and facets will be added.
*
* size_t DD_HIGHWATER
* size_t DD_LOWWATER
*    when a block has DD_WIGHWATER more bytes allocated, it is
*    reduced back to DD_LOWWATER
*/

/** maximal dimension we are willing to handle **/
#define MAXIMAL_ALLOWED_DIMENSION	100

/** initial number of vertices/facets. **/
#define DD_INITIAL_VERTEXNO	4092	// their numbmer fluctuates widely
#define DD_INITIAL_FACETNO	124	// their number increases by one

/** controlling too large memory blocks, in bytes **/
#define DD_HIGHWATER	((size_t)50e+6)   // 50M
#define DD_LOWWATER	((size_t)10e+6)   // 10M

/** asking space for 128 vertices and 4096 facets **/
#ifdef BITMAP_32		/* 32 bit bitmap blocks */
#define DD_VERTEX_ADDBLOCK	128
#define DD_FACET_ADDBLOCK	4
#else				/* 64 bit bitmaps blocks */
#define DD_VERTEX_ADDBLOCK	64
#define DD_FACET_ADDBLOCK	2
#endif

#ifdef USETHREADS
/************************************************************************
* Threads
*
* int create_threads(void)
*    initialize all threads. Return non-zero in case of an error.
*
* void stop_threads(void)
*    clear up threads
*/
int create_threads(void);
void stop_threads(void);

#endif /* USETHREADS */

/************************************************************************
* Statistics
*
* Several values are collected during the run of the DD algorithm,
*   which will be printed out when requested. Error conditions, such
*   as inconsistency, or out of memory, are also indicated here.
*
* void get_dd_vertexno(void)
*   computes the number of verticess of the most recent approximation as
*   well as how many of them is known to be a facet of the final
*   polyhedron (living_vertex_no and final_vertex_no). Keeping these 
*   values up-to-date at each iteration would be too expensive.
*/
typedef struct {
/** statistics **/
int iterations;		    /* number of iterations which is the same as
			       the number of facets added */
int facetno;                /* facets generated so far */
int vertexenquiries;	    /* number of times a new vertex was requested */
int probefacet;		    /* number of time facets were probed */
int vertices_allocated_no;  /* number of times vertex space was extended */
int vertices_allocated;	    /* total number of vertices allocated */
int facets_allocated_no;    /* number of times facet space was extended */
int facets_allocated;	    /* total number of facets allocated */
int vertex_compressed_no;   /* times vertex compression is called */
int vertex_pos;             /* last number of positive vertices */
int vertex_zero;            /* vertices adjacent to the current facet */
int vertex_neg;             /* negative vertices to be dropped */
int vertex_new;             /* total number of vertices added */
int max_vertices;	    /* maximal number of intermediate vertices */
int max_vertexadded;	    /* maximal number of vertices added in an iteration */
double avg_vertexadded;	    /* average number of vertices added */
double max_tests;	    /* maximal number of edge tests by iterations */
double avg_tests;	    /* average number of edge tests by iterations */
int living_vertex_no;	    /* filled by get_dd_vertexno() */
int final_vertex_no;	    /* filled by get_dd_vertexno() */
int memory_allocated_no;    /* number of times memory expanded */
size_t total_memory;	    /* total memory (in bytes) actually allocated */
size_t max_memory;          /* maximum memory allocated so far */
/** warning **/
int instability_warning;    /* number of warnings when recalculating facet eqs */
/** error conditions **/
int numerical_error;	    /* numerical error, data is inconsistent */
int out_of_memory;	    /* out of memory, cannot continue */
int data_is_consistent;     /* in case of memory shortage, indicate data consistency */
} DD_STATS;

extern DD_STATS dd_stats;
void get_dd_vertexno(void);

/************************************************************************
* Iterations of the double description algorithm
*
* int init_dd_structure(vertexno,facetno)
*    initialize DD algorithm structure allocating space to accomodate
*    the given number of vertices and facets, but at last DD_INITIAL_
*    many. Return value:
*     0:  initialization is succesful
*     1:  either dimension is too large, or out of memory.
*
* void init_dd(double v[0:dim-1])
*    initialize the DD algorithm, supplying the very first vertex.
*    Call init_dd_structure(0,0) before this routine.
*    The positive endpoints of coordinate axes (the *ideal* vertices) 
*    are assumed to be feasible solutions. Arguments:
*      v[0:dim-1]: coordinates of the first external vertex.
*
* int add_initial_vertex(int final,double coords[0..dim-1])
* int add_initial_facet(double coords[0..dim])
*    Load vertices and facets from an earlier partial computation.
*    First all facets, followed by the vertices. For a vertex 'final' 
*    indicates whether the vertex is final. Adjacency lists are
*    generated when adding a vertex.
*    Return value:
*     0:  OK
*     1:  some error, error message issued.
*
* int get_next_vertex(int from,double *v[0:dim-1])
*    Return the index and its coordinates of a vertex of the actual
*    approximation which is not marked as finial. When 'from' is non-
*    negative, start searching there. If from==-1, then return the 
*    smallest index, or a random index depending on the parameter
*    RandomVertex. Return -1 if no vertex was found.
*    
* void mark_vertex_as_final(vno)
*    Mark vertex with index vno as a vertex of the final polytope.
*
* void add_new_facet(double v[0:dim])
*    Specify the vertex which enlarges the inner approximation. When
*    the routine returns, check error conditions in dd_stats.
*/

/** initialize data structures with the first vertex **/
int init_dd_structure(int vertexno, int facetno);
void init_dd(const double *coords);

/** add initial vertex and facet **/
int add_initial_vertex(int final,const double coords[]);
int add_initial_facet(const double coords[]);

/** release vertex and facet adjacency lists **/
void free_adjacency_lists(void);

/** get next living but not final vertex **/
int get_next_vertex(int from, double *to);

/** mark the facet as final **/
void mark_vertex_as_final(int vno);
/** add a new facet **/
void add_new_facet(double *coords);

/************************************************************************
* Retrieving data, checking, recalculating
*
* int get_vertexnum()
* int get_facetnum()
*    Return the number of non-ideal vertices and facets.
*
* int probe_facet(double v[0:dim])
*    Return the facet score; facet with the highest score will be passed
*    to add_new_facet(). The number of vertices thrown away seems to be
*    a good heuristic.
*
* void recalculate_vertices(void)
*    Recalculate all vertices from the list of facets adjacent to it. 
*     When the routine returns, error conditions should be checked.
*/

/** actual number of vertices and facets */
int get_vertexnum(void); int get_facetnum(void); 

/** facet score, the higher the better **/
int probe_facet(double *coords);

/** recalculate verteices **/
void recalculate_vertices(void);

/************************************************************************
* Report the facets and vertices of the solution
*
* void print_vertex(report_type channel,int vno)
*    Report the coordinates of vertex 'vno' on the given channel. Use
*    fractional format if VertexAsFraction is set, otherwise a floating
*    number. Coordinates are separated by a space. A terminating
*    newline is added at the end.
*
* void print_vertices(report_type channel)
*    Report all vertices on the given channel. The coordinates of each
*    vertex is printed on a separate line. The first character is 'V'
*    or 'v'; coordinates are printed using print_vertex().
*
* void print_facet(report_type channel,double coodes[])
*    Report the facet equation. The coefficients are scaled making
*    them close to integers, if possible
*
* void print_facets(report_type channel)
*    Report all living facets of on the given channel using print_facet()
*
* void report_memory_usage(report_type ch, int force, char *prompt)
*    report memory usage of the inner approximation algorithm.
*
* void make_checkpoint(void)
*    create the next checkpoint file. Do not report any problem.
*
* void make_dump(void)
*    create a dump of vertices and facets when requested by a signal.
*
* int check_bitmap_consistency(void)
*    check number of adjacency numbers in the bitmaps. Return the
*    number of errors. (0 = no error)
*/

/** reporting vertices **/
void print_vertex(report_type channel, int vno);
void print_vertices(report_type channel); 

/** reporting facets **/
void print_facet(report_type channel,const double coords[]);
void print_facets(report_type channel);

/** report memory usage **/
void report_memory_usage(report_type channel, int force, const char *prompt);

/* create checkpoint files **/
void make_checkpoint(void);

/* create dump */
void make_dump(void);

/* bitmap consistency */
int check_bitmap_consistency(void);
int check_consistency(void);

/* EOF */


