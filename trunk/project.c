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

#define OUTPUT_INTERVAL 200	// Number of timesteps between each output to filex
#define MAX_COUNT 10000	        // Maximum number of iterations
#define TIMING_INTERVAL 5000     // Number of timesteps between timing stats
#define DISPLAY_INTERVAL 5     // Number of timesteps between each output to screen 
#define ANIMATION_INTERVAL 5   // Number of timesteps between each output for animation

/********** Variable Declarations **********/

MPI_Datatype MPI_STAR;

/********** Function Headers **********/

void initialize();
void do_gravitation();
void perform_calculations();

/********** Function Declarations **********/

/**
 * Initializes MPI environment.  Determines comm size and rank.  Provides
 * the loop for each timestep.  Prints to screen and to file as specified in
 * OUTOUT_INTERVAL and PRINTOUT_INTERVAL.
 */
int main(int argc, char* argv[])
{
	int count = 0;
	char str[100], cstr[100];
	unsigned long long start, end, s, e;

	// initialize MPI environment
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		printf("MPI intialization failed.\n");
		exit(1);
	}

	start = rdtsc();
	
	// get number of processors, current rank
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &my_size);
	
	// get filename from command line args
	if (argc < 3) {
		printf("USAGE: mpirun -np n [inputStarsFile] [darkMatterFile]\n");
		exit(1);
	} else { // get stars from file
	  initialize(argv[1], argv[2]);
	}
 
	while (count < MAX_COUNT) {
	        s = rdtsc();
		// sync up all processors
		if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
			printf("MPI_Barrier error in processor %d \n", my_rank);
			exit(1);
		}
		e = rdtsc();
		barrier_time += ((double)(e-s)) / CLOCK_RATE_KRATOS;

		perform_calculations();
		count++;

		// print out state if needed
		if (count % OUTPUT_INTERVAL == 0) {
		  strcpy(str, "./states");
		  sprintf(cstr, "%02d", my_size);
		  strcat(str, cstr);
		  strcat(str, "/outFile");
		  sprintf(cstr, "%04d", count);
		  strcat(str, cstr);
		  strcat(str, ".txt");
		  printStarInfo(str);
		}

		// print out positions for animations
		if (count % ANIMATION_INTERVAL == 0) {
		  strcpy(str, "./animations");
		  sprintf(cstr, "%02d", my_size);
		  strcat(str, cstr);
		  strcat(str, "/outFile");
		  sprintf(cstr, "%04d", count);
		  strcat(str, cstr);
		  strcat(str, ".txt");
		  printAnimations(str);
		}

		// print out step to screen
		if (my_rank == 0 && count % DISPLAY_INTERVAL == 0) 
		        printf("%d: Finished step %d.\n", my_rank, count);
		
		if (count % TIMING_INTERVAL == 0) 
			printTimingResults();
	}

	end = rdtsc();

	printf("%d: Total Time = %lf\n", my_rank, ((double)(end-start))/CLOCK_RATE_KRATOS);

	MPI_Finalize();	

	return 0;
}

/**
 * Call helper functions to load the data, create the array,
 * and otherwise set up the program's initial conditions.
 */
void initialize(char* fileName, char* darkMatter) 
{
  int lengths[3] = {1, 4, 1};
  MPI_Aint disp[3] = {0, 0, 16};
  MPI_Datatype types[3] = {MPI_LB, MPI_DOUBLE, MPI_UB};
	
  MPI_Type_create_struct(3, lengths, disp, types, &MPI_STAR);
  MPI_Type_commit(&MPI_STAR);
  
  // populate galaxy
  getStarInfo(fileName, darkMatter);	

}

/**
 * Applies gravitation to stars found in the "stars" array.  We sent the index of the 
 * star to the apply_gravitation function.
 */
void do_gravitation()
{
	int i;

	// apply gravitation due to stars and dark matter
	for (i=0;i<num_stars+num_dark;i++)
	        apply_gravitation(i);
}

/**
 * Performs all the calculations needed in each timestep.  This includes all
 * message passing (which passes around the galaxy arrays).
 */
void perform_calculations()
{
	int i;
	int round = 0, r_flag, s_flag, dest;
	MPI_Request recv_req, send_req;
	MPI_Status recv_status, send_status;
	unsigned long long r_start, r_end, s_start, s_end, start, end;

	// initialize acceleration
	for (i=0;i<num_stars;i++)
	{
		stars[i].x_acc = 0.0;
		stars[i].y_acc = 0.0;
		stars[i].z_acc = 0.0;
	}
	
	while (round < my_size) {
	  if (round != 0) { // if not the first round, wait for previous recieve to be done
	    do {
	      MPI_Test(&recv_req, &r_flag, &recv_status);
	    } while (r_flag == 0);
	    r_end = rdtsc();
	    message_time += ((double)(r_end-r_start)) / CLOCK_RATE_KRATOS;

	    galaxy = recv_array;
	  }
	  
	  // recieve array from the previous processor
	  if (round != my_size - 1) {
	    dest = (my_rank - 1 + my_size) % my_size;
	    MPI_Irecv(recv_array, num_stars + num_dark, MPI_STAR, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_req);
	    r_start = rdtsc();
	  }
	  // apply gravitation to the portion of current galaxy
	  start = rdtsc();
	  do_gravitation();
	  end = rdtsc();
	  computation_time += ((double)(end-start)) / CLOCK_RATE_KRATOS;
		
	   
	  if (round > 0) {
	    // need to make sure last send finished before sending again
	    do {
	      MPI_Test(&send_req, &s_flag, &send_status);
	    } while (s_flag == 0);    
	    s_end = rdtsc();
	    message_time += ((double)(s_end-s_start)) / CLOCK_RATE_KRATOS;
	  }
	  
	  if (round != my_size - 1) { // if not the last loop	    
	    // send block to next processor
	    dest = (my_rank + 1) % my_size;
	    MPI_Isend(galaxy, num_stars + num_dark, MPI_STAR, dest, 0, MPI_COMM_WORLD, &send_req);
	    s_start = rdtsc();
	    
	  }
	  round++;
	}

	// update velocities and accelerations
	for (i=0;i<num_stars;i++)
	  {
		stars[i].x_pos = stars[i].x_pos + stars[i].x_v + 0.5 * stars[i].x_acc;
		stars[i].y_pos = stars[i].y_pos + stars[i].y_v + 0.5 * stars[i].y_acc;
		stars[i].z_pos = stars[i].z_pos + stars[i].z_v + 0.5 * stars[i].z_acc;

		stars[i].x_v = stars[i].x_v + stars[i].x_acc;
		stars[i].y_v = stars[i].y_v + stars[i].y_acc;
		stars[i].z_v = stars[i].z_v + stars[i].z_acc;		
	  }

	// copy stars back into galaxy
	for (i=0; i < num_stars+num_dark;i++) {
	  galaxy[i].x_pos = stars[i].x_pos;
	  galaxy[i].y_pos = stars[i].y_pos;
	  galaxy[i].z_pos = stars[i].z_pos;
	  galaxy[i].mass = stars[i].mass;
	}

}

#endif
