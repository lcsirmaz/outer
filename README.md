# outer - a MultiObjective Linear Program (MOLP) solver

A MOLP is a linear program with more than one objective functions. The linear
*constraints* are arranged in a matrix form with c columns and r rows.
Columns correspond to variables x<sub>1</sub>, x<sub>2</sub>, . . . ,
x<sub>c</sub>, which are subject to the restrictions

<table><tbody><tr>
<td>L<sub>i</sub> &le; x<sub>i</sub> &le; U<sub>i</sub></td>
<td>where 1 &le; i &le; c</td>
</tr></tbody></table>

The lower bound L<sub>i</sub> can be -&#x221e; and the upper bound can be +&#x221e;.

For each row 1 &le; j &le; r the j-th constraint is

<table><tbody><tr>
<td>b<sub>j</sub> = a<sub>1,j</sub> x<sub>1</sub> + a<sub>2,j</sub> x<sub>2</sub> + . . . + a<sub>c,j</sub> x<sub>c</sub></td>
<td>l<sub>j</sub> &le; b<sub>j</sub> &le; u<sub>j</sub></td>
<td>where 1 &le; j &le; r</td>
</tr></tbody></table>

A *feasible solution* is a tuple **x** = &lt;x<sub>1</sub>, x<sub>2</sub>, . . . , x<sub>c</sub>&gt; 
which satisfies all constraints.

There are m &ge; 1 *objectives*, which are given as

<table><tbody><tr>
<td>y<sub>k</sub> = o<sub>1,k</sub> x<sub>1</sub> + o<sub>2,k</sub> x<sub>2</sub> + . . . + o<sub>c,k</sub> x<sub>c</sub></td>
<td>where 1 &le; k &le; m.</td>
</tr><tbody></table>

The m-dimensional objective vector **y** = &lt;y<sub>1</sub>, . . ., y<sub>m</sub>&gt; is 
*achievable* if there is a feasible solution **x** which yields exactly these objective values.

The *solution* of a MOLP is the list of all extremal objective vectors. For a 
minimizing (maximizing) MOLP the objective vector **y** is (weakly) *extremal* if there 
is no other achievable objective vector **y**' which would be coordinatewise &le; **y** 
(or &ge; when maximizing).

The task of MOLP solver is to find the collection of all extremal objective vectors.
These extremal vectors form the boundary of an m-dimensional convex polyhedral cone,
thus it suffices to find its *vertices*. This MOLP solver finds these vertices by iteratively 
computing outer approximations of the solution cone. In each iteration the 
approximating polytope is cut by a new facet (bounding hyperplane) of the
final cone. The time required for an iteration varies widely from a 
couple of microseconds to several days. After each iteration the solver
checks if the process has been interrupted by a SIGUSR1 signal.
If it has, it switches to a quick and dirty method which
might find further extremal vectors (but not necessarily all of them).


#### RESTRICTIONS

The linear constraints must have a feasible solution. Moreover, in a
minimizing MOLP all objectives must be bounded from below, and in a
maximizing MOLP all objectives must be bounded from above.


#### ALGORITHM

The algorithm is an implementation of Benson's outer approximation mathod,
see Benson, H.P.: *An Outer Approximation Algorithm for Generating All Efficient Extreme Points in the Outcome Set of a Multiple Objective Linear Programming Problem*. 
Journal of Global Optimization 13, 1–24 (1998). [https://doi.org/10.1023/A:1008215702611](https://doi.org/10.1023/A:1008215702611).

#### USAGE

The program is invoked as

    outer [options] <vlp-file>

or for the threaded version,

    outerth [options] <vlp-file>

The only obligatory argument is the file name which contains the description
of the problem in vlp format. Accepted options are

| Option | meaning |
|:-------|:--------|
| `-h`          | display a short help and quit |
| `--help`      | display all options |
| `--help=<topic>`  | help on one of the following topics: input, output, config, <br> exit, signal, checkpoint, resume, boot, vlp |
| `--help=output`  | describe the output format |
| `--version`   | version and copyright information |
| `--dump`      | dump the default config file and quit |
| `--config=<config-file>` or <br> `-c <config-file>`  | read configuration from the given file <br> use `--dump` to show the default config file |
| `-o <file>`  | save results (both vertices and facets) to \<file\> |
| `-ov <file>` | save vertices to \<file\> |
| `-of <file>` | save facets to \<file\> |
| `-oc <stub>` | file stub for checkpoint files |
| `--name=NAME` or <br> `-n NAME`    | specify the problem name |
| `--boot=<vertex-list>` | start the algorithm with these facets |
| `--resume=<chk-file>` | resume computation from a checkpoint file |
| `-m[0..3]`   | set message level: 0: none, 1: errors, 2: all, 3: verbose |
| `-q`         | quiet, same as `-m0`. Implies `--PrintStatistics=0` |
| `-p T`       | progress report in every T seconds (default: T=5) |
| `-p 0`       | no progress report |
| `-y+`        | report extremal solutions (vertices) immediately when generated (default) |
| `-y-`        | do not report extremal solutions when generated |
| `--KEYWORD=value` | change value of a config keyword |

#### CONFIGURATION PARAMETERS

Fine tuning the algorithm and the underlying scalar LP solver, and 
specifying the amount and type of information saved is done by giving
values of several keywords. Each keyword has a default value, which is
overwritten by the values in the config file (if specified), and those
values are overwritten by the `--KEYWORD=value` command line options.
Change tolerances with great care.

|**Algorithm parameters**| |
|:--------|:------------|
|`RandomVertex=1`<br>&nbsp; | 0 = no, 1 = yes <br>  pick the next vertex to be passed to the oracle randomly. |
|`RandomIdealPoint=1` <br>&nbsp; | 0 = no, 1 = yes <br> choose the ideal direction randomly rather than using (1,1,1,...,1). |
|`ExactVertexEq=0`<br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br>  when a vertex is created, recompute the coordinates immediately from the set of its adjacent facets. |
|`RecalculateVertices=100`<br>&nbsp;<br>&nbsp; | non-negative integer <br> after that many iterations recalculate all vertices from the set of its adjacent facets. The number should be zero (meaning never), or at least 5. |
|`CheckConsistency=0`<br>&nbsp;<br>&nbsp; | non-negative integer <br> after this many iterations check the consistency of the data structure against numerical errors. The number should be zero (meaning never), or at least 5. |
|`ExtractAfterBreak=1`<br>&nbsp;<br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> when the program receives a `SIGUSR1` signal or reaches the memory limit, continue extracting new vertices by asking the oracle about every current vertex of the approximating cone. Second signal aborts this post-processing. |
|`MemoryLimit=0`<br>&nbsp;<br>&nbsp; | non-negative integer <br> upper limit for memory allocation in Mbytes. When reaching this limit, stop processing as if received a `SIGUSR1` signal. Zero means no limit, otherwise it must be at least 100. |
|`TimeLimit=0`<br>&nbsp;<br>&nbsp; | non-negative integer <br> upper limit for running time in seconds. When reaching this limit, stop processing as if received a `SIGUSR1` signal. Zero means no limit, otherwise it must be at least 60 seconds. |
|`FacetPoolSize=0`<br>&nbsp;<br>&nbsp;<br>&nbsp; | non-negative integer <br> size of the facet pool: use the facet which discards the largest number of existing vertices. Should be zero (don't use it), or at least 5. Using facet pool adds more oracle calls, but can simplify the approximating polytope. |
|`CheckPoint=10000`<br>&nbsp;<br>&nbsp;<br>&nbsp;<br>&nbsp; | positive integer <br> frequency (in seconds) for creating checkpoint files when the option `-oc <filestub>` is given. The filename is got from the stub by appending `NNN.chk` where `NNN` starts with 000 and increases. The computation can be resumed from a checkpoint file by calling outer with `--resume=<checkpoint>` The value should be more than 500. |
|`OracleCallLimit=1`<br>&nbsp;<br>&nbsp; | non-negative integer <br> the maximal number of unsuccessful oracle calls during an iteration when filling the facet pool. Zero means no limit; otherwise should be less than 100. |
|`Threads=0`<br>&nbsp;<br>&nbsp; | non-negative integer <br> number of threads to use; should be less than 64. Zero means use as many threads as are available; 1 means don't use threads. |
|**Oracle parameters**| |
|`OracleMessage=1`<br>&nbsp; | 0 = quiet, 1 = error, 2 = on, 3 = verbose <br> oracle (glpk) message level. |
|`OracleMethod=0`<br>&nbsp;  | 0 = primal, 1 = dual <br> the LP method used by the oracle. |
|`OraclePricing=1`<br>&nbsp; | 0 = standard, 1 = steepest edge <br> the LP pricing method. |
|`OracleRatioTest=1`<br>&nbsp; | 0 = standard, 1 = Harris' two pass <br> the LP ratio test. |
|`OracleTimeLimit=20` <br>&nbsp; | non-negative integer <br> time limit for each oracle call in seconds; 0 means unlimited. |
|`OracleItLimit=10000` <br>&nbsp; | non-negative integer <br> iteration limit for each oracle call; 0 means unlimited. |
|`OracleScale=1` <br>&nbsp; | 0 = no, 1 = yes <br> scale the constraint matrix; helps numerical stability. |
|`ShuffleMatrix=1` <br>&nbsp; | 0 = no, 1 = yes <br> shuffle rows and columns of the constraint matrix randomly. |
|`RoundFacets=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> when the oracle reports a result facet, round its coordinates to the nearest rational number with small denominator. |
|**Reporting**| |
|`MessageLevel=2` <br>&nbsp;<br>&nbsp; | 0 = quiet, 1 = error, 2 = on, 3 = verbose <br> report level; quiet means no messages at all. The command line option `-m[0..3]` overrides this value. |
|`Progressreport=5` <br>&nbsp;<br>&nbsp; | non-negative integer <br> minimum time between two progress reports in seconds. Should be zero for no progress reports, or at least 5. The command line option `-p T` overrides this value. |
|`VertexReport=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> report each extremal vertex immediately when it is found. The command line option `-y-` (no) or `-y+` (yes) overrides the value defined here. |
|`FacetReport=0` <br>&nbsp; | 0 = no, 1 = yes <br> report each facet immediately when generated. |
|`MemoryReport=0` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = at the end, 2 = always <br> report the size and location, whenever they change, of the memory blocks storing the combinatorial data structure. |
|`VertexAsFraction=1` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = yes <br> if possible, print (and save) vertex coordinates as fractions with small denominators rather than floating point numerals. |
|`PrintStatistics=1` <br>&nbsp; | 0 = no, 1 = yes <br> print resources used (number of iterations, edge tests, etc.) when the program stops. |
|`PrintParams=1` <br>&nbsp; | 0 = no, 1 = yes <br> print algorithm parameters which are not equal to their default values. |
|`PrintVertices=2` <br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> print (again) all known vertices when the program terminates. |
|`PrintFacets=0` <br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> print all known facets when the program terminates. |
|`SaveVertices=2` <br>&nbsp;<br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> when the program terminates, save known vertices to the file specified after command line option `-o`. For file specified after `-ov` both 0 and 1 means &quot;save on normal exit only&quot;. |
|`SaveFacets=1` <br>&nbsp;<br>&nbsp;<br>&nbsp; | 0 = no, 1 = on normal exit only, 2 = always <br> when the program terminates, save known facets to the file specified after <br> the command line option `-o`. For the file specified after `-of` both 0 and 1 means &quot;save on normal exit only&quot;. |
|**Tolerances**|**Change these values with care**|
|`PolytopeEps=1.3e-8` <br>&nbsp; | positive real number <br> a facet and a vertex are considered adjacent if their distance is smaller than this value. |
|`ScaleEps=3e-9` <br>&nbsp;<br>&nbsp; | positive real number <br> coefficients in the scaled facet equation are rounded to the nearest integer if they are closer to it than this value. |
|`LineqEps=8e-8` <br>&nbsp;<br>&nbsp; | positive real number <br> when solving a system of linear equations for vertex coordinates, a coefficient smaller than this is considered to be zero. |
|`RoundEps=1e-9` <br>&nbsp;<br>&nbsp; | positive real number <br> if facets coordinates reported by the oracle are to be rounded (`RoundFacets=1`), this is the tolerance in the rounding algorithm. |
|`VertexRecalcEps=1e-6` <br>&nbsp;<br>&nbsp; | positive real number <br> when recalculating vertices, report numerical instability if the new and old coordinates differ at least this much. |

#### EXIT STATUS

When the program terminates, the exit status is zero when the problem has
been solved successfully; otherwise it indicates the failure condition:

|**Exit value**| | 
|:------|:-----|
|0    | problem solved |
|1    | error before starting the algorithm, such as missing argument, wrong or missing data file, out of memory, etc. |
|2    | there is no feasible solution for the MOLP problem |
|3    | the solution space is unbounded in some objective direction |
|4    | error while executing the algorithm (oracle failure, numerical error, etc) |
|5    | program execution was interrupted by a `SIGUSR1` signal |


#### SIGNALS

When the program receives a `SIGUSR1` signal, it stops the main cycle of
iterations, and switches to a &quot;quick and dirty&quot; method to generate
additional extremal solutions. (Actually, the oracle is asked about all
vertices of the actual approximation to check whether it is extremal.)
The method might miss extremal solutions, so the result is not known 
(can not) be complete. A second `SIGUSR1` signal aborts this post-processing.

When receiving a `SIGUSR2` signal, the program creates a snapshot file containing the 
vertices and facets of the actual approximation. Similarly to checkpoint
files, the snapshot file can be used to resume the computation from that point.

The name of the snapshot file is created as follows. If a checkpoint stub is provided
after the `-oc` option, then `000.dmp` is appended to the stub. If no
`-oc` option is provided but there is a `-o` output file, then the 
extension of that filename is replaced by `dmp`. If neither `-oc` nor `-o`
option is present, then no snapshot is created. When receiving a second
`SIGUSR2` signal, the same filename is used: the previous content is silently
overwritten.

Signals are checked only when an iteration step has been completed, that is a
new facet has been added to the approximation. (These are the points when
the data is guaranteed to be consistent.) Depending on the complexity of the
problem and the number of intermediate vertices and facets, the time
between these points vary from milliseconds to several days.

#### USING MULTIPLE CORES

The threaded version of `outer` can use all available cores on the computer 
to speed up the computation. Multiple cores can help in the combinatorial part, 
when the work is distributed approximately uniformly among them. The LP
part (oracle calls) is inherently single-threaded and cannot be speed up
this way. Thus LP intensive problems (large number of columns and rows
and only a few objectives) will not gain too much from multiple threads.
Specifying more cores than available comes with a significant performance
penalty. 


#### COMPILATION

The program uses glpk, the GNU Linear Program Kit, for solving scalar LP problems.
Changing to the directory `OUTER/src`, the following command compiles **outer** 
linking the glpk routines:

    gcc -O3 -W -o outer *.c -lm -lglpk

The threaded version requires defining `USETHREADS` and adding the flag `-pthread`:

    gcc -O3 -W -I -o outerth -DUSETHREADS -pthread *.c -lm -lglpk

#### AUTHOR

Laszlo Csirmaz, <csirmaz@ceu.edu>


