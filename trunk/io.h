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
#include "timing.h"

/********** Variable Definitions **********/

/********** Variable Declarations **********/
int io_count;             
int b_count;
int comp_count;
int msg_count;
double barrier_time = 0.0;      // Keeps track of all the time spent in barriers
double computation_time = 0.0;  // Keeps track of all the time spent in computation
double message_time = 0.0;      // Keeps track of all the time spent in message passing
double max_x;
double max_y;
double max_z;

/********** Function Headers **********/

void getStarInfo(char* fileName, char* darkfile);
int printStarInfo(char* fileName, int print_velocity);
void printTimingResults();

/********** Function Declarations **********/

/**
 * Parses through input file to retrieve all star info.  All information is
 * put into an array of stars - stars and galaxy.
 *
 * stars holds the stars tied to this processor
 * galaxy holds the current stars (at the start, stars and galaxy will be the same)
 *
 * INPUTS: char* file name
 * OUTPUT: none
 *
 * Times are printed to timing_io_results.txt.
 */
void getStarInfo(char* fileName, char* darkMatter) {
  MPI_File file, darkfile;
  int i, offset;
  int line_length = 133; // 18 chars per field, 7 fields, 6 separating commas, 1 newline
  double x_p, y_p, z_p, x_vl, y_vl, z_vl, m;
  char buf[1024], str[100], cstr[100];
  MPI_Status status;
  unsigned long long start, end;

  // reset counters
  io_count = b_count = comp_count = msg_count = 0;


  // get size of dark matter file
  num_dark = 17071;
  if (num_dark % my_size != 0)
      num_dark -= num_dark % my_size;
  num_dark /= my_size;

  start = rdtsc();
  
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

  num_stars = NUMBER_OF_STARS / my_size;

  if (my_rank == 0)
    printf("%d: We have %ld stars, each processor has %d stars and %d dark matter 'stars'.\n",  my_rank, NUMBER_OF_STARS, num_stars, num_dark);
 
  galaxy = (star*) malloc((num_stars + num_dark) * sizeof(star));
  recv_array = (star*) malloc((num_stars + num_dark) * sizeof(star));
  stars = (star*) malloc((num_stars + num_dark) * sizeof(star));

  for (i = 0; i < num_stars; i++) {
    offset = my_rank * line_length * num_stars + i * line_length + 13;

    MPI_File_read_at(file, offset, buf, line_length, MPI_CHAR, &status);

    sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &x_p, &y_p, &z_p, &x_vl, &y_vl, &z_vl, &m);

    // malloc all galaxy arrays
    (stars[i]).x_pos = (galaxy[i]).x_pos = x_p;
    (stars[i]).y_pos = (galaxy[i]).y_pos = y_p;
    (stars[i]).z_pos = (galaxy[i]).z_pos = z_p;
    (stars[i]).x_v   = (galaxy[i]).x_v   = x_vl;
    (stars[i]).y_v   = (galaxy[i]).y_v   = y_vl;
    (stars[i]).z_v   = (galaxy[i]).z_v   = z_vl;
    (stars[i]).mass  = (galaxy[i]).mass  = m; 
  }
  
  // close file
  MPI_File_close(&file);
 
  end = rdtsc();

  // open dark matter file
  if (MPI_File_open(MPI_COMM_WORLD, darkMatter, MPI_MODE_RDONLY, MPI_INFO_NULL, &darkfile) != MPI_SUCCESS) {
    printf("%d: Error in opening file.\n", my_rank);
    exit(1);
  }

  line_length = 58;
  for (i = 0; i < num_dark; i++) {
    offset = my_rank * line_length * num_dark + i * line_length;

    MPI_File_read_at(file, offset, buf, line_length, MPI_CHAR, &status);

    sscanf(buf, "%lf,%lf,%lf\n", &x_p, &y_p, &z_p);

    // malloc all galaxy arrays
    (stars[i + num_stars]).x_pos = (galaxy[i + num_stars]).x_pos = x_p;
    (stars[i + num_stars]).y_pos = (galaxy[i + num_stars]).y_pos = y_p;
    (stars[i + num_stars]).z_pos = (galaxy[i + num_stars]).z_pos = z_p;
    (stars[i + num_stars]).x_v   = (galaxy[i + num_stars]).x_v   = 0;
    (stars[i + num_stars]).y_v   = (galaxy[i + num_stars]).y_v   = 0;
    (stars[i + num_stars]).z_v   = (galaxy[i + num_stars]).z_v   = 0;
    (stars[i + num_stars]).mass  = (galaxy[i + num_stars]).mass  = 4; 
  }
  
  // close file
  MPI_File_close(&darkfile);

  strcpy("./timing");
  sprintf(cstr, "%02d", my_size);
  strcat(str, cstr);
  strcat(str, "/timing_io_results.txt");

  // write out timing results
  if (MPI_File_open(MPI_COMM_WORLD, str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    printf("%d: Error in opening file.\n", my_rank);
    exit(1);
  }

  offset = io_count * my_size * 20 + my_rank * 20;
  
  sprintf(buf, "%04d,%014.9lf\n", my_rank, ((double)(end-start))/CLOCK_RATE_BGL);
 
  MPI_File_write_at(file, offset, buf, 20, MPI_CHAR, &status);
  
  // close file
  MPI_File_close(&file);
  
  io_count++;
}

/**
 * Goes through an array of stars and outputs their current state into
 * a file.  Times are printed to timing_results.txt.
 *
 * INPUTS: char* file name
 * OUTPUT: 1 = error in creating file
 *         0 = success
 */
int printStarInfo(char* fileName, int print_velocity) {

  // get file, open it
  MPI_File file;
  char buf[1024], str[100], cstr[100]; 
  int i, offset, line_length;
  MPI_Status status;
  unsigned long long start, end;

  if (MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    fprintf(stderr, "Cannot create %s.\n", fileName);
    exit(1);
  }

  if (my_rank == 0)
    printf("%d: Printing out to file: %s.\n", my_rank, fileName);

  if (print_velocity == 0)
    line_length = 57;
  else
    line_length = 133;

  start = rdtsc();
  // print out galaxy size
  if (my_rank == 0) {
    sprintf(buf, "%012ld\n", NUMBER_OF_STARS);      
    MPI_File_write_at(file, 0, buf, 13, MPI_CHAR, &status);
  }

  // loop through array and output state of each star
  for (i = 0; i < num_stars; i++) {
    offset = my_rank * line_length * num_stars + i * line_length + 13;
    
    if (print_velocity == 1)
    sprintf(buf, "% 018lf,% 018lf,% 018lf,% 018lf,% 018lf,% 018lf,% 018lf\n", stars[i].x_pos, stars[i].y_pos, stars[i].z_pos, stars[i].x_v, stars[i].y_v, stars[i].z_v, stars[i].mass);
    else
      sprintf(buf, "% 018lf,% 018lf,% 018lf\n", stars[i].x_pos, stars[i].y_pos, stars[i].z_pos);

    MPI_File_write_at(file, offset, buf, line_length, MPI_CHAR, &status);
  }
  end = rdtsc();

  // close file 
  MPI_File_close(&file);

  strcpy("./timing");
  sprintf(cstr, "%02d", my_size);
  strcat(str, cstr);
  strcat(str, "/timing_io_results.txt");
  
  // write out timing results
  if (MPI_File_open(MPI_COMM_WORLD, str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    printf("%d: Error in opening file.\n", my_rank);
    exit(1);
  }

  offset = io_count * my_size * 20 + my_rank * 20;
  
  sprintf(buf, "%04d,%014.9lf\n", my_rank, ((double)(end-start))/CLOCK_RATE_BGL);

  MPI_File_write_at(file, offset, buf, 20, MPI_CHAR, &status);
  
  // close file
  MPI_File_close(&file);

  io_count++;
  
  return 0;
}
/**
 * Prints barrier, computation, and message results to file (per processor)
 */
void printTimingResults() 
{
  int offset;
  char buf[100];
  MPI_File file;
  MPI_Status status;

  strcpy("./timing");
  sprintf(cstr, "%02d", my_size);
  strcat(str, cstr);
  strcat(str, "/barrier_results.txt");

  // BARRIER
  if (MPI_File_open(MPI_COMM_WORLD, str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    printf("%d: Error in opening barrier file.\n", my_rank);
    exit(1);
  }

  offset = b_count * my_size * 20 + my_rank * 20;
  
  sprintf(buf, "%04d,%014.9lf\n", my_rank, barrier_time);

  MPI_File_write_at(file, offset, buf, 20, MPI_CHAR, &status);
  
  // close file
  MPI_File_close(&file);

  b_count++;

  // COMPUTATION
  if (MPI_File_open(MPI_COMM_WORLD, "./timing/computation_results.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    printf("%d: Error in opening computation file.\n", my_rank);
    exit(1);
  }

  offset = b_count * my_size * 20 + my_rank * 20;
  
  sprintf(buf, "%04d,%014.9lf\n", my_rank, computation_time);

  MPI_File_write_at(file, offset, buf, 20, MPI_CHAR, &status);
  
  // close file
  MPI_File_close(&file);

  comp_count++;

  strcpy("./timing");
  sprintf(cstr, "%02d", my_size);
  strcat(str, cstr);
  strcat(str, "/message_results.txt");

  // MESSAGE
  if (MPI_File_open(MPI_COMM_WORLD, str, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
    printf("%d: Error in opening file.\n", my_rank);
    exit(1);
  }

  offset = b_count * my_size * 20 + my_rank * 20;
  
  sprintf(buf, "%04d,%014.9lf\n", my_rank, message_time);

  MPI_File_write_at(file, offset, buf, 20, MPI_CHAR, &status);
  
  // close file
  MPI_File_close(&file);

  msg_count++;


}
#endif
