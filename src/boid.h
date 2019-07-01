/*
 *  boid.h
 *  BoidTest
 *
 *  Created by maruyama on 11/04/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MASSIVESWARM_BOID_H
#define MASSIVESWARM_BOID_H

#include "vector3D.h"

class Boid {
public:
	Boid();
	Boid(Vector3D _p);
	Boid(Vector3D _p, Vector3D _v);
	void set(double x, double y, double z, double vx, double vy, double vz);
	void get(double* x, double* y, double* z, double* vx, double* vy, double* vz);
	//void set_serialized_data(double* buffer);
	void set_serialized_data(float* buffer);
	//void get_serialized_data(double* buffer);
	void get_serialized_data(float* buffer);
	bool isInsideArea(Boid &_b, double _sighe_distance, double _sight_angle);
    //bool isInsideSeparationArea(Boid &_b);
    //bool isInsideAlignmentArea(Boid &_b);
	//bool isInsideCohesionArea(Boid &_b);

	unsigned id;
	Vector3D position;
	Vector3D velocity;

	int grid_index[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
};

typedef struct _boid {
    unsigned int id;
    float x, y, z;
    float vx, vy, vz;
} boid;

#endif
