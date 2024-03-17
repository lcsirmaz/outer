/** glp_oracle.c  --  facet separation oracle **/

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
* This part of OUTER realizes the facet separation oracle using glpk
* (Gnu Linear Program Kit). The solution to the request is in a random
* direction or in the direction (1,1,1,...,1).
* OracleData structure is used to communicate with the calling routine.
*/

#include "glp_oracle.h"
#include "report.h"
#include "params.h"

Oracle_t OracleData;

/**********************************************************************
* The problem dimensions, question and answer space
*   int vcols   number of columns in the constraint matrix
*   int vrows   number of rows in the constraint matrix
*   int vobjs   dimension of the objective space: number of objectives
*   double vvertex[vobjs]
*               the question direction
*   double vfacet[vobjs]
*               the answer vertex which minimizes the facet direction
*/

#define vcols	PARAMS(ProblemColumns)
#define vrows	PARAMS(ProblemRows)
#define vobjs	PARAMS(ProblemObjects)
#define vvertex OracleData.overtex
#define vfacet	OracleData.ofacet

/*==================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 	// sqrt
#include "glpk.h"

/**********************************************************************
* The glpk LP object and parameter structure
*   glp_prob *P    the glpk problem object
*   glp_smcp parm  the glpk parameter structure
*/

static glp_prob *P=NULL;
static glp_smcp parm; /* glp parameter structure */

/**********************************************************************
* Read the constraint matrix and objectives from a vlp file
*
* int load_vlp()
*   open and read the problem from the vlp file. Return value:
*     0: file read, memory allocated, dimensions stored,
*        glpk LP object P is initialized and ready to use.
*     1: some error; errors are reported as R_fatal. The vlp file
*        is not necessarily closed; the program should abort.
*
*/

/*---------------------------------------------------------------------
* void perm_array(int len, int array[1..len])
*   make a random permutation of the array
*/
static inline int mrandom(int v)
{ return v<=1?0 : (random()&0x3fffffff)%v; }
static inline void perm_array(int len, int arr[/* 1:len */])
{int i,j,t;
    for(i=1;i<len;i++){
        j=i+mrandom(len+1-i);
        t=arr[i];arr[i]=arr[j];arr[j]=t;
    }
}
/*---------------------------------------------------------------------
* Storage for the constraint matrix, the objectives and the shuffle
*   arrays. The matrix and the shuffle arrays will be freed after
*   loading the problem into 'P'.
*
* double M(row,col)             temporary storage for the constraint
*                               matrix; indices go from 1
* int vlp_rowidx[1:rows+objs]   shuffling rows, temporary
* int vlp_colidx[1:cols+objs+1] shuffling columns, temporary
* int vlp_objidx[1:objs]        objective indices in rows
* int vlp_slackobj[1:objs]      slack variables in columns
* int lambdaidx                 lambda column index
* double vlp_lambda[1:objs]     obj coeffs in the lambda column
*
* int allocate_vlp(row,cols,objs)
*   allocate memory for vfacet, vvertex, and temporal storage.
*   Return value: -1: out of memory (should not happen)
*/

static double *vlp_M;		/* temporary storage for M */
static int *vlp_rowidx;		/* temporary row permutation */
static int *vlp_colidx;		/* temporary column permutation */
static int *vlp_objidx;		/* object indices */
static int *vlp_slackobj;	/* slack variables */
static double *vlp_lambda;	/* lambda column */
static int lambda_idx;		/* lambda column index */

/* M has size row+objs times cols+objs+1 */
#define M(r,c)		vlp_M[(r)+((c)-1)*(vrows+vobjs)-1]

#define xalloc(ptr,type,size)	\
    (ptr=(type*)calloc(size,sizeof(type)))==NULL

static int allocate_vlp(int rows, int cols, int objs)
{ int i;
    vrows=rows; vcols=cols; vobjs=objs; // store these values in PARAMS()
    if(xalloc(vfacet,double,objs+2) ||
       xalloc(vvertex,double,objs+1) ||
       xalloc(vlp_objidx,int,objs+1) ||
       xalloc(vlp_lambda,double,objs+1) ||
       xalloc(vlp_slackobj,int,objs+1) ||
       xalloc(vlp_M,double,(rows+objs)*(cols+objs+1)) ||
       xalloc(vlp_rowidx,int,rows+objs+1) ||
       xalloc(vlp_colidx,int,cols+objs+2))
         return -1; // out of memory
    // indirect inidices to rows and columns
    for(i=0;i<=rows+objs;i++) vlp_rowidx[i]=i;
    for(i=0;i<=cols+objs+1;i++) vlp_colidx[i]=i;
    if(PARAMS(ShuffleMatrix)){ // permute randomly
        perm_array(rows+objs,vlp_rowidx);
        perm_array(cols+objs+1,vlp_colidx);
    }
    // objective rows and slack columns, and the lambda column
    for(i=1;i<=objs;i++) vlp_objidx[i]=vlp_rowidx[i+rows];
    for(i=1;i<=objs;i++) vlp_slackobj[i]=vlp_colidx[i+cols];
    lambda_idx=vlp_colidx[cols+objs+1];
    return 0;
}

/*---------------------------------------------------------------------
* Character input from the vlp file
* int MAX_LINELEN = 80
*   maximum line length expected in a vlp file
* char inpline[MAX_LINELEN]
*   the next line read form the vlp file
*
* int nextline(FILE *f)
*   read the next line to inpline[] from the open stream 'f'. Ignore
*       leading spaces and merge spaces. Convert A-Z to a-z.
*   Return value:
*        1: next line is in inpline[]
*        0: EOF
*
* int vlp_type_ok(char ctrl,int parno)
*   checks the ctrl char against the supplied number of args:
*    'f' (free) parno==0
*    'u','l','s' (upper,lower, fixed) parno==1
*    'd' (both) parno==2
*
* int glp_type(char ctrl)
*   returns the glpk version of the ctrl char, or -1 if error.
*
* double ideal_direction()
*    either +1, or a random number between 0 and 1 with non-zero
*    low-precision bits to make the computation more stable
*/

/* read a single line from a VLP file */
#define MAX_LINELEN	80 /* maximal length of a file */

static char inpline[MAX_LINELEN+1]; /* contains the next vlp line */

/* read the next vlp line to inpline */
static int nextline(FILE *f)
{int i,sp,ch;
    i=0;sp=0; memset(inpline,0,MAX_LINELEN+1);
    while((ch=getc(f))>=0){
        if(ch=='\n'){ if(i==0){sp=0; continue;} return 1; }
        if(ch==' '||ch=='\t'){sp=1; continue;}
        if(ch<=0x20 || ch>126) continue; /* ignore these characters */
        if(sp && i>0){if (i<MAX_LINELEN){inpline[i]=' '; i++;} }
        if('A'<=ch && ch<='Z') ch +='a'-'A'; /* upper case => lower case */
        sp=0; if(i<MAX_LINELEN){inpline[i]=ch; i++; }
    }
    /** EOF **/
    return i>0?1:0;
}

/* direction character and the number of following real numbers */
static int vlp_type_ok(char ctrl, int parno)
{   switch(ctrl){
  case 'f': return parno==2;
  case 'u': case 'l': case 's': return parno==3;
  case 'd': return parno==4;
    }
    return 0;
}

/* glpk version of directions */
static int glp_type(char ctrl)
{   switch(ctrl){
  case 'f': return GLP_FR;
  case 'u': return GLP_UP;
  case 'l': return GLP_LO;
  case 's': return GLP_FX;
  case 'd': return GLP_DB;
    }
    return -1;
}

/* either 1.0 or a random number with non-zero least significant bits */
static double ideal_direction(void)
{   return
      PARAMS(RandomIdealPoint)
       ? sqrt(0.699+0.3*((double)random())/((double)0x7fffffff))
       : 1.0;
}

/* read a vlp problem from a file as an LP instance */
int load_vlp(void)
{FILE *f; int rows,cols,objs; int i,j,cnt; double p,b1,b2; char ctrl;
 double dir=1.0;
    f=fopen(PARAMS(VlpFile),"r");
    if(!f){
        report(R_fatal,"Cannot open vlp file %s for reading\n",PARAMS(VlpFile));
        return 1;
    }
    rows=cols=0;
    P = glp_create_prob();
    while(nextline(f)) switch(inpline[0]){ /* read next input line */
       case 'c': // comment line, print out nonempty comment lines before the first p line
                 if(rows==0 && inpline[1]){
                     report(R_warn,"C%s\n",inpline+1);
                 }
                 continue;
       case 'e': break; // end
       case 'p': if(rows>0){
                    report(R_fatal,"read_vlp: second p line in %s:\n   %s\n",
                              PARAMS(VlpFile),inpline); return 1; 
                 }
                 dir=+1.0; PARAMS(Direction)=0;
                                    /**   <rows> <cols> <> <objs> <> **/
                 cnt=sscanf(inpline,"p vlp min %d %d %*d %d %*d",&rows,&cols,&objs);
                 if(cnt==0){
                    dir=-1.0; PARAMS(Direction)=1;
                    cnt=sscanf(inpline,"p vlp max %d %d %*d %d %*d",&rows,&cols,&objs);
                 }
                 if(cnt!=3 || rows<=1 || cols<=1 || objs<1 ){
                    report(R_fatal,"read_vlp: wrong p line in %s\n   %s\n",
                               PARAMS(VlpFile),inpline); return 1;
                 }
                 if(allocate_vlp(rows,cols,objs)){
                     report(R_fatal,"read_vlp: out of memory for %s:\n   %s\n",
                               PARAMS(VlpFile),inpline); return 1;
                 }
		 glp_add_cols(P,cols+objs+1); glp_add_rows(P,rows+objs);
                 glp_set_col_bnds(P,lambda_idx,GLP_FR,0.0,0.0); // lambda is free
                 continue;
       case 'j': // j <col> [ f || l <val> | u <val> | d <val1> <val2> | s <val> ]
                 if(rows==0){
                    report(R_fatal,"read_vlp: j line before p in %s\n  %s\n",
                               PARAMS(VlpFile),inpline); return  1;
                 }
                 b1=b2=0.0;
                 cnt=sscanf(inpline,"j %d %c %lg %lg",&j,&ctrl,&b1,&b2);
                 if(cnt<2 || cols<j || j<1 || !vlp_type_ok(ctrl,cnt)){
                    report(R_fatal,"read_vlp: wrong j line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 if(cnt<4) b2=b1; // GLP_UP uses b2 as the value
                 glp_set_col_bnds(P,vlp_colidx[j],glp_type(ctrl),b1,b2);
                 continue;
       case 'i': if(rows==0){
                    report(R_fatal,"read_vlp: i line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 // i <row> [ f | l <val> | u <val> | d <val1> <val2> | s <val> ]
                 b1=b2=0.0;
                 cnt=sscanf(inpline,"i %d %c %lg %lg",&i,&ctrl,&b1,&b2);
                 if(cnt<2 || rows<i || i<1 || !vlp_type_ok(ctrl,cnt)){
                    report(R_fatal,"read_vlp: wrong i line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 if(cnt<4) b2=b1; // GLP_UP uses b2 as the value
                 glp_set_row_bnds(P,vlp_rowidx[i],glp_type(ctrl),b1,b2);
                 continue;
       case 'a': if(rows==0){
                    report(R_fatal,"read_vlp: a line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 cnt=sscanf(inpline,"a %d %d %lg",&i,&j,&p);
                 if(cnt!=3 || rows<i || i<1 || cols<j || j<1){
                    report(R_fatal,"read_vlp: wrong a line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 M(vlp_rowidx[i],vlp_colidx[j])=p; // store it
                 continue;
       case 'o': if(rows==0){
                    report(R_fatal,"read_vlp: o line before p in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 cnt=sscanf(inpline,"o %d %d %lg",&i,&j,&p);
                 if(cnt!=3 || objs<i || i<1 || cols<j || j<1){
                    report(R_fatal,"read_vlp: wrong o line in %s\n   %s\n",
                                PARAMS(VlpFile),inpline); return 1;
                 }
                 M(vlp_objidx[i],vlp_colidx[j])=dir*p; // store it
                 continue;
       default: report(R_fatal,"read_vlp: unknown line in %s\n  %s\n",
                           PARAMS(VlpFile),inpline); return 1;
    }
    fclose(f);
    if(rows==0){
       report(R_fatal,"read_vlp: no 'p' line in %s\n",PARAMS(VlpFile)); return 1; 
    }
    /* the vlp file has been read; set the glpk LP instance */
    // slack variables; they are fixed to be zero, to be changed later to GLP_LO
    for(j=1;j<=objs;j++) glp_set_col_bnds(P,vlp_slackobj[j],GLP_FX,0.0,0.0);
    // diagonal entries of the slack variables are 1.0
    for(j=1;j<=objs;j++) M(vlp_objidx[j],vlp_slackobj[j])=1.0;
    // finally upload constraints into P
    for(i=0;i<=rows+objs;i++) vlp_rowidx[i]=i; // index file
    for(j=1;j<=cols+objs+1;j++){
        if(j!=lambda_idx) glp_set_mat_col(P,j,rows+objs,vlp_rowidx,&M(1,j)-1);
    }
    free(vlp_colidx); free(vlp_rowidx); free(vlp_M);
    // negtive lambda coefficients are in vlp_lambda
    for(j=1;j<=objs;j++) vlp_lambda[j]= -1.0*ideal_direction();
    glp_set_mat_col(P,lambda_idx,objs,vlp_objidx,vlp_lambda);
    for(j=1;j<=objs;j++) vlp_lambda[j] *= -1.0; // revert to positive values
    // LP objective: minimize lambda
    glp_set_obj_coef(P,lambda_idx,1.0);
    glp_set_obj_dir(P,GLP_MIN);
    // preprocess the constraint matrix
    if(PARAMS(OracleMessage)<2) glp_term_out(GLP_OFF);
    glp_sort_matrix(P);
    if(PARAMS(OracleScale)) glp_scale_prob(P,GLP_SF_AUTO);
    glp_adv_basis(P,0);  // make this optimization
    glp_term_out(GLP_ON); // enable messages form glpk
    return 0;
}

/**********************************************************************
* void set_oracle_parameters(void)
*   Set the LP solver parameters from the configuration:
*    verbosity: 0: no, 1: error; 2: normal, 3: verbose
*    output frequency: indicate that the LP solver is working
*    method: primal, dual
*    pricing: standard or steepest edge
*    ratio test: standard of Harris
*    iteration limit
*    time limit (in seconds)
*/

static void set_oracle_parameters(void)
{   glp_init_smcp(&parm);
    // verbosity
    switch(PARAMS(OracleMessage)){
      case  0: parm.msg_lev=GLP_MSG_OFF; break;
      case  1: parm.msg_lev=GLP_MSG_ERR; break;
      case  2: parm.msg_lev=GLP_MSG_ON; break;
      default: parm.msg_lev=GLP_MSG_ALL; break;
    }
    // method, pricing, ratio test
    parm.meth = GLP_PRIMAL;	// PRIMAL,DUAL,DUALP
    if(PARAMS(OracleMethod)) parm.meth=GLP_DUAL;
    parm.pricing = GLP_PT_STD;
    if(PARAMS(OraclePricing)) parm.pricing=GLP_PT_PSE;
    parm.r_test = GLP_RT_STD;	// HAR
    if(PARAMS(OracleRatioTest)) parm.r_test=GLP_RT_HAR;
    // iteration and time limit
    parm.it_lim = 100000;	// iteration limit
    if(PARAMS(OracleItLimit)>=1000) parm.it_lim=PARAMS(OracleItLimit);
    if(PARAMS(OracleItLimit)==0) parm.it_lim=0; // no limit
    parm.tm_lim = 10000;	// time limit in milliseconds
    if(PARAMS(OracleTimeLimit)>=5) parm.tm_lim=1000*PARAMS(OracleTimeLimit);
    if(PARAMS(OracleTimeLimit)==0) parm.tm_lim=0; // no limit
}

/*********************************************************************&
* char *glp_status_msg(int code)
* char *glp_return_msg(int code)
*  return the verbatim error message corresponding to the glpk
*  code (to be printed out).
*/

static char *glp_status_msg(int stat)
{static char *statmsg[] = {
"the problem is undefined",		// GLP_UNDEF
"solution is feasible",			// GLP_FEAS
"solution is infeasible",		// GLP_INFEAS
"the problem has no feasible solution",	// GLP_NOFEAS
"solution is optimal",			// GLP_OPT
"the problem is unbounded",		// GLP_UNBND
};
    if(1<=stat && stat<=6) return statmsg[stat-1];
    return "unknown solution status";
}

static char *glp_return_msg(int retval)
{static char *retmsg[] = {
"invalid basis",			// GLP_EBADB	 *
"singular matrix",			// GLP_ESING	 *
"ill-conditioned matrix",		// GLP_ECOND	 *
"invalid bounds",			// GLP_EBOUND
"solver failed",			// GLP_EFAIL	 *
"objective lower limit reached",	// GLP_EOBJLL
"objective upper limit reached",	// GLP_EOBJUL
"iteration limit exceeded",		// GLP_EITLIM	 *
"time limit exceeded",			// GLP_ETMLIM	 *
"no primal feasible solution",		// GLP_ENOPFS
"no dual feasible solution",		// GLP_ENODFS
"root LP optimum not provided",		// GLP_EROOT
"search terminated by application",	// GLP_ESTOP
"relative mip gap tolerance reached",	// GLP_EMIPGAP
"no primal/dual feasible solution",	// GLP_ENOFEAS
"no convergence",			// GLP_ENOCVG
"numerical instability",		// GLP_EINSTAB
"invalid data",				// GLP_EDATA
"result out of range",			// GLP_ERANGE
};
    if(1<=retval && retval<=0x13) return retmsg[retval-1];
    return "unknown error";
}

/**********************************************************************
* int initialize_oracle(void)
*  Check if the consistency of the loaded vlp prolem; add the ideal
*  points if the are missing
*
* init get_initial_vertex(void)
*   Create an ouside point to OracleData.overtex for the first approximation.
*
* int ask_oracle(void)
*  The question and the answer is provided in OracleData. Return values:
*/

#include "round.h" /* round_to */

/* count the number and measure time of LP calls */
static int oracle_calls=0;
static unsigned long oracle_time=0ul;

/* call the glpk simplex solver twice if necessary */
#include <sys/time.h> 
static int call_glp(void)
{int ret; struct timeval tv; unsigned long starttime;
    oracle_calls++;
    if(gettimeofday(&tv,NULL)==0){
        starttime=tv.tv_sec*1000 + (tv.tv_usec+500u)/1000u;
    } else {starttime=0ul;}
    ret=glp_simplex(P,&parm);
    if(ret==GLP_EFAIL){ // give it a second chance
        if(PARAMS(OracleMessage)<2) glp_term_out(GLP_OFF);
        glp_adv_basis(P,0);
        glp_term_out(GLP_ON);
        oracle_calls++;
        ret=glp_simplex(P,&parm);
    }
    if(gettimeofday(&tv,NULL)==0)
        oracle_time += (tv.tv_sec*1000 + (tv.tv_usec+500u)/1000u)-starttime;
    return ret;
}

int initialize_oracle(void)
{int i,j,ret;
    set_oracle_parameters();
    // compute the maximum along the i-th objective
    glp_set_obj_dir(P,GLP_MAX);
    for(i=1;i<=vobjs;i++){
        // set the objective rows bouunds
        for(j=1;j<=vobjs;j++)
            glp_set_row_bnds(P,vlp_objidx[j],i==j ? GLP_FX : GLP_FR,0.0,0.0);
        ret=call_glp();
        if(ret){
            report(R_fatal,"The oracle says: %s (obj %d)\n",glp_return_msg(ret),i);
            return ORACLE_FAIL;
        }
        ret=glp_get_status(P);
        if(ret==GLP_OPT){ // bounded in this direction
            glp_set_col_bnds(P,vlp_slackobj[i],GLP_LO,0.0,0.0);
        } else if(ret!=GLP_UNBND){
            report(R_fatal,"The oracle says: %s (%d)\n",glp_status_msg(ret),ret);
            return ret==GLP_NOFEAS ? ORACLE_EMPTY : ORACLE_FAIL;
        } // otherwise OK
    }
    // from now on we compute minimum
    glp_set_obj_dir(P,GLP_MIN);
    return ORACLE_OK;
}

int get_initial_vertex(void)  // create an outside vertex
{int i,j,ret; double d;
    for(i=1;i<=vobjs;i++){
       // set objective rows bounds all free except for a single fixed
       for(j=1;j<=vobjs;j++) 
            glp_set_row_bnds(P,vlp_objidx[j],i==j ? GLP_FX : GLP_FR,0.0,0.0);
       ret=call_glp();
       if(ret){
            report(R_fatal,"The oracle says: %s (obj %d)\n",glp_return_msg(ret),i);
            return ORACLE_FAIL;
        }
        ret=glp_get_status(P);
        if(ret != GLP_OPT){
            report(R_fatal,"The oracle says: %s (obj %d)\n",glp_status_msg(ret),i);
            return ret==GLP_UNBND ? ORACLE_UNBND : ORACLE_FAIL;
        }
        // recover the solution
        d=glp_get_obj_val(P)*vlp_lambda[i];
        if(PARAMS(RoundFacets))round_to(&d);
        vvertex[i-1]=d;
    }
    return ORACLE_OK;
}

/* int ask_oracle() 
*   ask oracle about OracleData.overtex, return OracleData.ofacet as the 
*   separating supporting hyperplane
*/
int ask_oracle(void)
{int i,ret; double lambda,d;
    // set the objective bound to the values in the vertex
    for(i=1;i<=vobjs;i++)
        glp_set_row_bnds(P,vlp_objidx[i],GLP_UP,0.0,vvertex[i-1]);
    ret=call_glp();
    if(ret){
        report(R_fatal,"The oracle says: %s (%d)\n",glp_return_msg(ret),ret);
        // one can continue if  ret==GLP_EITLIM || ret==GLP_ETMLIM
        return ORACLE_FAIL;
    }
    ret=glp_get_status(P);
    if(ret == GLP_UNBND) // inside
        return ORACLE_UNBND;
    if(ret != GLP_OPT){
        report(R_fatal,"The oracle says: %s (%d)\n",glp_status_msg(ret),ret);
        return ORACLE_FAIL; 
    }
    lambda=glp_get_obj_val(P); // optimal value of lambda
    /* if lambda <=0 the vertex is inside the polytope
     * the boundary point: << vvertex[i-1]+lambda*vlp_lambda[i] >>
     * the facet equation is the negative of the dual solution */
    if(lambda < PARAMS(PolytopeEps)) return ORACLE_UNBND;
    for(i=1;i<=vobjs;i++) vfacet[i-1]=-glp_get_row_dual(P,vlp_objidx[i]);
    // normalize the equation to sum up to 1.0
    d=0.0;
    for(i=0;i<vobjs;i++) d += vfacet[i];
    for(i=0;i<vobjs;i++) vfacet[i] /= d;
    for(i=0;i<vobjs;i++){
        d=vfacet[i];
        if(PARAMS(RoundFacets)) round_to(&d);
        // all coefficients must be >=0
        if(d<0.0){
            if(d<-PARAMS(PolytopeEps)){
               report(R_fatal,"Numerical instability, facet coeff[%d]=%g<0\n",i+1,d);
               return ORACLE_FAIL;
            }
            d=0.0;
        }
        vfacet[i]=d; }
    // the optimal solution is on the supporting hyperplane
    d=0.0;
    for(i=1;i<=vobjs;i++){
        d+=vfacet[i-1]*(vvertex[i-1]+lambda*vlp_lambda[i]);
    }
    if(PARAMS(RoundFacets)) round_to(&d); 
    vfacet[vobjs]=-d;
    return ORACLE_OK;
}

/**********************************************************************
* Get oracle statistics
*
* void get_oracle_stat(int *no, init *it, unsigned long *time)
*   return the number of LP calls and time. Use the 
*   undocumented glpk call to get the number of iterations.
*/

void get_oracle_stat(int *no, int *it, unsigned long *time)
{   *no=oracle_calls; *it=glp_get_it_cnt(P); 
    *time=(oracle_time+5ul)/10ul; // in 0.01 seconds
}

/* EOF */

