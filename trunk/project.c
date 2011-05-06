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

#ifndef GALAXY_SIM_PROJECT_C
#define GALAXY_SIM_PROJECT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "star.h"
#include "io.h"
#include "timing.h"

/********** Variable Definitions **********/

#define OUTPUT_INTERVAL 100	// Number of timesteps between each 
                                // output of the current simulation state
#define MAX_COUNT 500	// maximum number of iterations

/********** Variable Declarations **********/

/********** Function Headers **********/

void initialize();
void update_galaxy();

/********** Function Declarations **********/

int main(int argc, char* argv[])
{
  int my_rank, my_size;
  int count = 0, i, num_stars;
  char str[100], cstr[100];

  // initialize MPI environment
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    printf("MPI intialization failed.\n");
    exit(1);
  }
  
  // get number of processors, current rank
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &my_size);
  
  // get filename from command line args
  if (my_rank == 0) {
    if (argc < 2) {
      printf("USAGE: mpirun -np _ [filename] [DEBUG] \n");
      exit(1);
    } else { // get stars from file
      initialize(argv[1]);
    }
  }

  while (count < MAX_COUNT) {   
    // sync up all processors
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
      printf("MPI_Barrier error in processor %d \n", my_rank);
      exit(1);
    }

    // calculate new positions, update star array
    num_stars = NUMBER_OF_STARS / my_size;

    for (i = 0; i < num_stars; i++)
      apply_gravitation(my_rank * num_stars + i);
    
    count++;

    // print out state if needed
    if ((my_rank == 0) && (count % OUTPUT_INTERVAL == 0)) {
      printf("%d: Printing out status at interval %d.\n", my_rank, count);
      strcpy(str, "./states/outFile");
      sprintf(cstr, "%d", count);
      strcat(str, cstr);
      strcat(str, ".txt");
      printStarInfo(str);
    }
  }
  
  MPI_Finalize();
  
  return 0;
}

/**
 * Call helper functions to load the data, create the array,
 * and otherwise set up the program's initial conditions.
 */
void initialize(char* fileName) 
{
  getStarInfo(fileName);
}

/**
 * Set the state of the galaxy (galaxy) to the new state after computation (new galaxy)
 */
void update_galaxy()
{
  star** old = galaxy;
  galaxy = new_galaxy;

  new_galaxy = old;

  memset(new_galaxy, 0, sizeof(old));
}

#endif
