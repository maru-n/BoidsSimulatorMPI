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
	void set_serialized_data(double* buffer);
	void get_serialized_data(double* buffer);
	bool isInsideArea(Boid &_b, double _sighe_distance, double _sight_angle);
    //bool isInsideSeparationArea(Boid &_b);
    //bool isInsideAlignmentArea(Boid &_b);
	//bool isInsideCohesionArea(Boid &_b);

	Vector3D position;
	Vector3D velocity;
};

#endif
