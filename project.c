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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "star.h"
#include "io.h"
#include "timing.h"

/********** Variable Definitions **********/

#define OUTPUT_INTERVAL 100  // Number of timesteps between each output 
                             // of the current simulation state

/********** Variable Declarations **********/

int my_rank;
int my_size;
MPI_Status status;
star* galaxy;


/********** Function Headers **********/

void initialize();
star get_star(int index);


/********** Function Declarations **********/

int main(int argc, char* argv[])
{
  // initialize MPI environment
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    printf("MPI intialization failed.\n");
    exit(1);
  }

  // get number of processors, current rank
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &my_size);

  initialize();






  MPI_Finalize();

  return 0;
}

/**
 * Call helper functions to load the data, create the array,
 * and otherwise set up the program's initial conditions.
 */
void initialize() 
{
  galaxy = getStarInfo("initialStarConfig.dat");
}

/**
 * Return the star at a given index in the array of stars
 */
star get_star(int index)
{
  return galaxy[index];
}