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

#ifndef GALAXY_SIM_STAR_H
#define GALAXY_SIM_STAR_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct
{
	//Components of star position
	double x_pos;
	double y_pos;
	double z_pos;
	//Components of star velocity
	double x_v;
	double y_v;
	double z_v;
	//Mass of star -- currently all stars are assumed to have the same mass
	double mass;
	//Acceleration components -- used only for gravitational calculations, not otherwise accessed directly
	double x_acc;
	double y_acc;
	double z_acc;
} star;

/********** Variable Definitions **********/

#define NUMBER_OF_STARS 10000      //The number of stars in the galaxy
#define STARS_IN_CLUSTER 1         //See description of cluster() function
#define UPDATE_INTERVAL 50         //Number of timesteps between updates of closest stars
#define GRAVITATION_DISTANCE 300   //Only stars within this maximum distance (in parsecs) will be taken into account for gravitation

/********** Variable Declarations **********/

double closest_cluster_stars[NUMBER_OF_STARS/STARS_IN_CLUSTER][STARS_IN_CLUSTER];  //Keeps track of star indices for faster cluster computation
star* galaxy[NUMBER_OF_STARS];
star* new_galaxy[NUMBER_OF_STARS];

/********** Function Headers **********/

star* cluster(star* cluster_stars[]);
double distance(int self_index, int other_index);
int* get_stars_within_range(int origin);
int* get_closest_stars(int origin);
void apply_gravitation(int origin);
star* force_of_gravity(int self_index, int other_index);
star* get_star(int index);
void set_star(int index, star* self);

/********** Function Declarations **********/

/**
 * Group the stars into clusters of CLUSTER_OF_STARS "virtual stars" to make calculations on the stars less memory-intensive when passing data around
 * via MPI.
 * 
 * INPUT: array of stars to be clustered
 * OUTPUT: star containing agreggate star cluster data
 */
star* cluster(star* cluster_stars[])
{
	star* collective;
	double sum_of_masses = 0.0;
	double sum_of_x = 0.0;
	double sum_of_y = 0.0;
	double sum_of_z = 0.0;
	int x;
	for(x=0;x<STARS_IN_CLUSTER;x++)
	{
		sum_of_masses += cluster_stars[x]->mass;
		sum_of_x += cluster_stars[x]->x_pos * cluster_stars[x]->mass;
		sum_of_y += cluster_stars[x]->y_pos * cluster_stars[x]->mass;
		sum_of_z += cluster_stars[x]->z_pos * cluster_stars[x]->mass;
	}
	sum_of_x /= sum_of_masses;
	sum_of_y /= sum_of_masses;
	sum_of_z /= sum_of_masses;
	collective->x_pos = sum_of_x;
	collective->y_pos = sum_of_y;
	collective->z_pos = sum_of_z;
	return collective;	
}

/**
 * Calculate the distance between two given stars
 * 
 * INPUT: the two stars under consideration
 * OUTPUT: 3-dimensional distance between the input stars as a double
 */
double distance(int self_index, int other_index)
{
	star* self = (star*)get_star(self_index);
	star* other = (star*)get_star(other_index);
	double dist = pow(other->x_pos - self->x_pos,2) + pow(other->y_pos - self->y_pos,2) + pow(other->z_pos - self->z_pos,2);
	return sqrt(dist);
}

/**
 * Determine which stars are within GRAVITATION_DISTANCE of the origin star
 * 
 * INPUT: the index of the origin star
 * OUTPUT: an array of the indices of the stars no farther than GRAVITATION_DISTANCE away
 */
int* get_stars_within_range(int origin)
{
	int size = NUMBER_OF_STARS / 10;
	int num_stars = 0;
	int* list = (int*)malloc(sizeof(int) * size);
	
	int x;
	for(x=0;x<NUMBER_OF_STARS;x++)
	{
		if(x != origin && distance(origin,x) < GRAVITATION_DISTANCE)
		{
			list[num_stars++] = x;
			if(num_stars >= size)  //Expand the array if necessary
			{
				size += 1;
				list = (int*)realloc(list, sizeof(int) * size);
			}
		}
	}
	return list;
}

/**
 * Determine the closest STARS_IN_CLUSTER stars to the origin star
 * 
 * INPUT: the index of the origin star
 * OUTPUT: an array of the indices of the STARS_IN_CLUSTER closest stars
 */
int* get_closest_stars(int origin)
{
	int* cluster = (int*)malloc(sizeof(int) * STARS_IN_CLUSTER);
	int* list_of_stars = get_stars_within_range(origin);
	int x,y;
	int index = 0;
	for(x=0;x<STARS_IN_CLUSTER;x++)
	{
		double min = 100000.0;
		for(y=0;x<sizeof(list_of_stars)/sizeof(int);y++)  //Loop through the list of stars within range, finding the closest
		{
			double dist = distance(origin,list_of_stars[y]);
			if(list_of_stars[y] != -1 && dist < min)
			{
				min = dist;
				index = y;
			}
		}
		cluster[x] = list_of_stars[y];
		list_of_stars[y] = -1;  //Once a closest star is found, set the index in the list to -1 so it isn't found again
	}
	return cluster;
}

/**
 * Apply the net gravitational force to the star
 * 
 * INPUT: the index of the star under consideration
 */
void apply_gravitation(int origin)
{
	int* other_stars = get_stars_within_range(origin);
	int x;
	star* self = (star*)get_star(origin);
	self->x_acc = self->y_acc = self->z_acc = 0.0;
	for(x=0;x<sizeof(other_stars)/sizeof(int);x++)
	{
		star* temp = force_of_gravity(origin,x);
		self->x_acc += temp->x_acc;
		self->y_acc += temp->y_acc;
		self->z_acc += temp->z_acc;
	}
	set_star(origin, self);
}

/**
 * Determine the magnitude of the gravitational force of one star on another
 * 
 * INPUT: the indices of the two stars under consideration
 * OUTPUT: a "dummy" star holding the resulting acceleration data
 */
star* force_of_gravity(int self_index, int other_index)
{
	double G = 4.49734287 * pow(10.0,-9);  //Gravitational constant, using units of parsecs * (solar mass units)^-1 * (parsecs/millennium)^2
	star* self = (star*)get_star(self_index);
	star* other = (star*)get_star(other_index);
	star* storage;
	double r = distance(self_index, other_index);
	double force = G * self->mass * other->mass;
	force /= pow(r,3);
	storage->x_acc = force * (other->x_pos - self->x_pos);
	storage->y_acc = force * (other->y_pos - self->y_pos);
	storage->z_acc = force * (other->z_pos - self->z_pos);
	return storage;
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
	new_galaxy[index] = self;
}

#endif
