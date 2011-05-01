/*
 * GalaxySim
 * Copyright (c) 2011 Cody Garges <gargec@rpi.edu>
 *                    James McAtamney <mcataj@cs.rpi.edu>
 *                    Amanda Olyha <olyhaa@rpi.edu>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef GALAXY_SIM_IO_H
#define GALAXY_SIM_IO_H

#include "star.h"
#include <stdio.h>
#include <stdlib.h>

/********** Variable Definitions **********/

/********** Variable Declarations **********/

/********** Function Headers **********/

void getStarInfo(char* fileName, star* galaxy[]);
int printStarInfo(char* fileName, star* galaxy[]);

/********** Function Declarations **********/


/*
 * Parses through input file to retrieve all star info.  All information is
 * put into an array of stars which is returned.
 *
 * INPUTS: char* file name
 *         int number of stars
 * OUTPUT: array of stars
 */
void getStarInfo(char* fileName, star* galaxy[]) {

  // get file, open it
  FILE *file;
  file = fopen(fileName, "r");
  if (file == NULL) {
    fprintf(stderr, "Cannot open %s to read in.\n", fileName);
    return;
  }

  // parse through file, adding to array
  int i = 0;
  float x_p, y_p, z_p, x_vl, y_vl, z_vl;

  while (fscanf(file, "%f, %f, %f, %f, %f, %f", &x_p, &y_p, &z_p, &x_vl, &y_vl, &z_vl) != EOF) {
    star* myStar;
    myStar->x_pos = x_p;
    myStar->y_pos = y_p;
    myStar->z_pos = z_p;
    myStar->x_v = x_vl;
    myStar->y_v = y_vl;
    myStar->z_v = z_vl;
    // myStar.mass = m;

    galaxy[i] = myStar;
    i++;
  }

  // close file
  fclose(file);
}



/*
 * Goes through an array of stars and outputs their current state into
 * a file.
 *
 * INPUTS: char* file name
 *         star* array of stars
 *         int number of stars
 * OUTPUT: 1 = error in creating file
 *         0 = success
 */
int printStarInfo(char* fileName, star* galaxy[]) {

  // get file, open it
  FILE *file;
  file = fopen(fileName, "w");
  if (file == NULL) {
    fprintf(stderr, "Cannot create %s to write out.\n", fileName);
    return 1;
  }

  int i;

  // loop through array and output state of each star
  for (i = 0; i < NUMBER_OF_STARS; i++) {
    fprintf(file, "%f, %f, %f, %f, %f, %f\n", (galaxy[i])->x_pos, (galaxy[i])->y_pos, (galaxy[i])->z_pos, (galaxy[i])->x_v, (galaxy[i])->y_v, (galaxy[i])->z_v);
  }

  // close file
  fclose(file);

  return 0;
}

#endif
