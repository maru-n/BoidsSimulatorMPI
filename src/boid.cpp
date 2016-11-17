/*
 *  Created by maruyama on 11/04/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 */

#include "boid.h"
#include <math.h>

static const double PI = 6*asin( 0.5 );

Boid::Boid(){
}

Boid::Boid(Vector3d _p){
	position = _p;
}

Boid::Boid(Vector3d _p, Vector3d _v){
	position = _p;
	velocity = _v;
}

bool Boid::isInsideArea(Boid &_b, double _sighe_distance, double _sight_angle) {
	if( ((*this).position - _b.position).norm() < _sighe_distance ){
        if (_sight_angle >= 1.0) {
            return true;
        }
		Vector3d sightVec = (*this).velocity;
        Vector3d relativePos = _b.position - (*this).position;
        double th = acos(sightVec.dot(relativePos) / (sightVec.norm() * relativePos.norm()));
		if( th <= PI*_sight_angle ){
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

