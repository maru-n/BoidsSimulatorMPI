/*
 *  boid.h
 *  BoidTest
 *
 *  Created by maruyama on 11/04/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "parameter.h"
#include <Eigen/Core>

using Eigen::Vector3d;

class Boid {
public:
	Boid();
	Boid(Vector3d _p);
	Boid(Vector3d _p, Vector3d _v);

    bool isInsideArea(Boid &_b, double _sighe_distance, double _sight_angle);
    bool isInsideSeparationArea(Boid &_b);
    bool isInsideAlignmentArea(Boid &_b);
	bool isInsideCohesionArea(Boid &_b);

	Vector3d position;
	Vector3d velocity;
};


