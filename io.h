/*
 * GalaxySim
 * Copyright (c) 2011 Cody Garges <gargec@rpi.edu>, James McAtamney <mcataj@cs.rpi.edu>, Amanda Olyha <olyhaa@rpi.edu>
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
#include "mpi.h"
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
  MPI_File file;
  if (MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    fprintf(stderr, "Cannot open %s to read in.\n", fileName);
    return;
  }

  // parse through file, adding to array
  float x_p, y_p, z_p, x_vl, y_vl, z_vl;
  char buf[1024]; 
  int i, offset, numStars = NUMBER_OF_STARS / getMySize(), line_length = 181; // 19 chars per field, 7 fields, 6 separating commas
  MPI_Status status;

  for (i = 0; i <numStars; i++) {
    offset = getMyRank() * line_length * numStars + i;
    MPI_File_read_at(file, offset, buf, line_length, MPI_CHAR, &status);

    sscanf(buf, "%f,%f,%f,%f,%f,%f", &x_p, &y_p, &z_p, &x_vl, &y_vl, &z_vl);

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
  MPI_File_close(&file);
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
  MPI_File file;
  if (MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    fprintf(stderr, "Cannot create %s.\n", fileName);
    return 1;
  }

  char buf[1024]; 
  int i, offset, numStars = NUMBER_OF_STARS / getMySize(), line_length = 181;
  MPI_Status status;

  // loop through array and output state of each star
  for (i = 0; i <numStars; i++) {
    offset = getMyRank() * line_length * numStars + i;
    
    sprintf(buf, "%f,%f,%f,%f,%f,%f\n", (galaxy[i])->x_pos, (galaxy[i])->y_pos, (galaxy[i])->z_pos, (galaxy[i])->x_v, (galaxy[i])->y_v, (galaxy[i])->z_v);
    
    MPI_File_write_at(file, offset, buf, line_length, MPI_CHAR, &status);
  }

  // close file 
  MPI_File_close(&file);
  
  return 0;
}

#endif
