/*
 *  Created by maruyama on 11/04/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 */

#include "boid.h"
#include "parameter.h"
#include <math.h>

Boid::Boid(){
	Vector3D p;
	Vector3D v;
	position = p;
	velocity = v;
}

Boid::Boid(Vector3D _p){
	Vector3D v;
	position = _p;
	velocity = v;
}

Boid::Boid(Vector3D _p, Vector3D _v){
	position = _p;
	velocity = _v;
}

bool Boid::isInsideArea(Boid &_b, double _sighe_distance, double _sight_angle) {
	if( ( (*this).position - _b.position ).norm() < _sighe_distance ){
        if (_sight_angle >= 1.0) {
            return true;
        }
		Vector3D sightVec = (*this).velocity;
        Vector3D relativePos = _b.position - (*this).position;
        double th = double(acos((sightVec * relativePos) / (sightVec.norm() * relativePos.norm())));
		if( th <= M_PI*_sight_angle ){
			return true;
        }else{
            return false;
        }
	}
	return false;
}

bool Boid::isInsideSeparationArea(Boid &_b) {
    return isInsideArea(_b, SIGHT_DISTANCE_SEPARATION, SIGHT_ANGLE_SEPARATION);
}

bool Boid::isInsideAlignmentArea(Boid &_b) {
    return isInsideArea(_b, SIGHT_DISTANCE_ALIGNMENT, SIGHT_ANGLE_ALIGNMENT);
}

bool Boid::isInsideCohesionArea(Boid &_b) {
    return isInsideArea(_b, SIGHT_DISTANCE_COHESION, SIGHT_ANGLE_COHESION);
}

