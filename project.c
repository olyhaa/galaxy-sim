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

#define OUTPUT_INTERVAL 100	// Number of timesteps between each output of the current simulation state
#define MAX_COUNT 500	// maximum number of iterations

/********** Variable Declarations **********/

double* recv_array = NULL;
int my_rank;
int my_size;
int num_stars;
star** stars;

/********** Function Headers **********/

void initialize();
void update_galaxy();

/********** Function Declarations **********/

int main(int argc, char* argv[])
{
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
	if (argc < 2) {
		printf("USAGE: mpirun -np _ [filename] \n");
		exit(1);
	} else { // get stars from file
		initialize(argv[1]);
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
	int count = 10;
	int lengths[3] = {1, 10, 1};
	MPI_Aint disp = {0, 0, 40};
	MPI_Datatype types[3] = {MPI_LB, MPI_DOUBLE, MPI_UB};
	MPI_Datatype MPI_STAR;
	MPI_Type_create_struct(3, lengths, disp, types, &MPI_STAR);
	MPI_Type_commit(&MPI_STAR);
	getStarInfo(fileName);
	num_stars = NUMBER_OF_STARS/my_size;
	recv_array = malloc(num_stars * sizeof(double));
	stars = malloc(num_stars * sizeof(double));
}

void do_gravitation(star* stars[], star* galaxy[])
{
	int i;
	star* temp;
	for (i=0;i<num_stars;i++)
	{
		temp = apply_gravitation(galaxy[i], galaxy);
		stars[i]->x_acc += temp->x_acc;
		stars[i]->y_acc += temp->y_acc;
		stars[i]->z_acc += temp->z_acc;
	}
}

void perform_calculations()
{
	int i,k;
	int tag = 0;
	MPI_Request req[2];
	for (i=0;i<num_stars;i++)
	{
		stars[i]->x_acc = 0.0;
		stars[i]->y_acc = 0.0;
		stars[i]->z_acc = 0.0;
	}
	//First, apply gravitation to our own portion of the galaxy
	do_gravitation(stars, galaxy);
	//Now send our part of the B matrix to processor my_rank+1, my_rank+2, ..., my_size-1, 0, 1, ..., my_rank-1 while receiving from the corresponding processors
	//Send/Recv is reversed by even/odd processor number to avoid deadlocks
	//Once a B matrix is received, do the necessary multiplication, then pass it along
	for(k=0;k<my_size-1;k++)
	{
		if(my_rank%2 == 0)
		{
			MPI_Isend(&galaxy[0],num_stars,MPI_STAR,(my_rank+1)&(my_size-1),tag,MPI_COMM_WORLD,&req[0]);
			MPI_Irecv(&recv_array[0],num_stars,MPI_STAR,(my_rank-1)&(my_size-1),tag,MPI_COMM_WORLD,&req[1]);
		}
		else
		{
			MPI_Irecv(&recv_array[0],num_stars,MPI_STAR,(my_rank-1)&(my_size-1),tag,MPI_COMM_WORLD,&req[1]);
			MPI_Isend(&galaxy[0],num_stars,MPI_STAR,(my_rank+1)&(my_size-1),tag,MPI_COMM_WORLD,&req[0]);
		}
		MPI_Waitall(2,req,MPI_STATUSES_IGNORE);
		//Swap galaxy and recv_array pointers, moving received data in recv_array to galaxy in preparation for calculations
		star** temp = galaxy;
		galaxy = recv_array;
		recv_array = temp;
		do_gravitation(stars, galaxy);
	}
}

#endif
