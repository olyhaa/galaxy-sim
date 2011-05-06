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

#include <stdio.h>
#include <stdlib.h>
#include "star.h"
#include "mpi.h"

/********** Variable Definitions **********/

/********** Variable Declarations **********/

/********** Function Headers **********/

void getStarInfo(char* fileName);
int printStarInfo(char* fileName);

/********** Function Declarations **********/

/**
 * Parses through input file to retrieve all star info.  All information is
 * put into an array of stars which is returned.
 *
 * INPUTS: char* file name
 * OUTPUT: none
 */
void getStarInfo(char* fileName) {
  MPI_File file;
  int my_rank, my_size, i, offset, numStars;
  int line_length = 133; // 18 chars per field, 7 fields, 6 separating commas, 1 newline
  double x_p, y_p, z_p, x_vl, y_vl, z_vl, m;
  char buf[1024];
  MPI_Status status;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &my_size);

  // get file, open it
  if (MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    printf("%d: Error in opening file.\n", my_rank);
    exit(1);
  }

  if (my_rank == 0)
    printf("%d: Opened file: %s.\n", my_rank, fileName);

  // get the number of stars
  MPI_File_read_at(file, 0, buf, 12, MPI_CHAR, &status);
  sscanf(buf, "%012ld", &NUMBER_OF_STARS);      
  
  if (NUMBER_OF_STARS % my_size != 0)
    NUMBER_OF_STARS -= NUMBER_OF_STARS % my_size;

  // barrier - to make sure NUMBER_OF_STARS has been set before continuing
  if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
    printf("MPI_Barrier error in processor %d \n", my_rank);
    exit(1);
  }

  numStars = NUMBER_OF_STARS / my_size;

  if (my_rank == 0)
    printf("%d: We have %ld stars, each processor has %d stars.\n",  my_rank, NUMBER_OF_STARS, numStars);
 
  galaxy = (star**) malloc(numStars * sizeof(star*));
  new_galaxy = (star**) malloc(numStars * sizeof(star*));

  for (i = 0; i < numStars; i++) {
    offset = my_rank * line_length * numStars + i * line_length + 13;

    MPI_File_read_at(file, offset, buf, line_length, MPI_CHAR, &status);

    sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &x_p, &y_p, &z_p, &x_vl, &y_vl, &z_vl, &m);

    galaxy[i] = (star*) malloc(sizeof(star));
    new_galaxy[i] = (star*) malloc(sizeof(star));
   
    (galaxy[i])->x_pos = x_p;
    (galaxy[i])->y_pos = y_p;
    (galaxy[i])->z_pos = z_p;
    (galaxy[i])->x_v = x_vl;
    (galaxy[i])->y_v = y_vl;
    (galaxy[i])->z_v = z_vl;
    (galaxy[i])->mass = m; 
  }
  
  // close file
  MPI_File_close(&file);
}

/**
 * Goes through an array of stars and outputs their current state into
 * a file.
 *
 * INPUTS: char* file name
 * OUTPUT: 1 = error in creating file
 *         0 = success
 */
int printStarInfo(char* fileName) {

  // get file, open it
  MPI_File file;
  char buf[1024]; 
  int my_rank, my_size, i, offset, numStars, line_length = 133;
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &my_size);

  if (MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    fprintf(stderr, "Cannot create %s.\n", fileName);
    exit(1);
  }

  numStars = NUMBER_OF_STARS / my_size;

  // print out galaxy size
  if (my_rank == 0) {
    sprintf(buf, "%012ld\n", NUMBER_OF_STARS);      
    MPI_File_write_at(file, 0, buf, 13, MPI_CHAR, &status);
  }

  // loop through array and output state of each star
  for (i = 0; i < numStars; i++) {
    offset = my_rank * line_length * numStars + i * line_length + 13;
    
    sprintf(buf, "% 018lf,% 018lf,% 018lf,% 018lf,% 018lf,% 018lf,% 018lf\n", (galaxy[i])->x_pos, (galaxy[i])->y_pos, (galaxy[i])->z_pos, (galaxy[i])->x_v, (galaxy[i])->y_v, (galaxy[i])->z_v, (galaxy[i])->mass);
    
    MPI_File_write_at(file, offset, buf, line_length, MPI_CHAR, &status);
  }

  // close file 
  MPI_File_close(&file);
  
  return 0;
}

#endif
