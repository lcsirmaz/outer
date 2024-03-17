/** outer.h  -- bounds and default values **/

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

/***********************************************************************
* the outerr approximation algorithm
*
* int outer(void)
*    execute the outer approximation algorithm for the prepared data.
*    Return value:
*      0   program terminated normally
*      1   data error before the algorithm started
*      2   problem unbounded
*      3   problem unfeasible
*      4   error during algorithm execution
*      5   interrupted, postprocessing terminated normally
*      6   interrupted, error in postprocessing
*      7   interrupted, postprocessing aborted
*/

int outer(void);

/* EOF */

