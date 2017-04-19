/*
 *  boid.h
 *  BoidTest
 *
 *  Created by maruyama on 11/04/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef parameter_H
#include "parameter.h"
#endif

//#include <glut/glut.h>
//#include <GLUI/GLUI.h>
//#include "/Users/ikeg/Downloads/glui-2.36/src/include/GL/glui.h"
//#include "vector3D.h"
#include "../vector3D.h"

class Boid {
public:
	Boid();
	Boid(Vector3D _p);
	Boid(Vector3D _p, Vector3D _v);

	//void draw();
    bool isInsideArea(Boid _b, double _sighe_distance, double _sight_angle);
    bool isInsideSeparationArea(Boid _b);
    bool isInsideAlignmentArea(Boid _b);
	bool isInsideCohesionArea(Boid _b);

	Vector3D position;
	Vector3D velocity;

	bool live;
};


