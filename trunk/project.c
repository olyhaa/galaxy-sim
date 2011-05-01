/*
 * GalaxySim
 * Copyright (c) 2011 Cody Garges <gargec@rpi.edu>
 *										James McAtamney <mcataj@cs.rpi.edu>
 *										Amanda Olyha <olyhaa@rpi.edu>
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
#include <string.h>
#include "star.h"
#include "io.h"
#include "timing.h"

/********** Variable Definitions **********/

#define OUTPUT_INTERVAL 100	// Number of timesteps between each output of the current simulation state
#define MAX_COUNT 50000000	 // maximum number of iterations

/********** Variable Declarations **********/

int my_rank;
int my_size;
MPI_Status status;
star* galaxy[NUMBER_OF_STARS];


/********** Function Headers **********/

void initialize();
star* get_star(int index);
void set_star(int index, star* self);

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

	// make sure the number of processors is valid
	if (NUMBER_OF_STARS % my_size != 0) {
		printf("The number of stars must be divisible by the processor world size.\n");
		exit(1);
	}

	// read in stars from file
	if (my_rank == 0)
		initialize();

	int count = 0;
	while (count < MAX_COUNT) {

		// sync up all processors
		if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
			printf("MPI_Barrier error in processor %d \n", my_rank);
			exit(1);
		}
		int i;
		// calculate new positions, update star array
		for (i = 0; i < NUMBER_OF_STARS / my_size; i++)
			apply_gravitation(my_rank * my_size + i);
	 
		count++;
		
		// print out state if needed
		if ((my_rank == 0) && (count % UPDATE_INTERVAL == 0)) {
		        char str[100], cstr[100];
			strcpy(str, "outFile");
			sprintf(cstr, "%d", count);
			strcat(str, cstr);
			strcat(str, ".dat");
			printStarInfo(str, galaxy);
		}
	}

	MPI_Finalize();

	return 0;
}

/**
 * Call helper functions to load the data, create the array,
 * and otherwise set up the program's initial conditions.
 */
void initialize() 
{
	getStarInfo("initialStarConfig.dat", galaxy);
}

/**
 * Return the star at a given index in the array of stars
 * 
 * INPUT: index of desired star
 */
star* get_star(int index)
{
	return galaxy[index];
}

/**
 * Set the star at index to new value
 * 
 * INPUT: index of desired star, desired star
 */
void set_star(int index, star* self)
{
	galaxy[index] = self;
}
