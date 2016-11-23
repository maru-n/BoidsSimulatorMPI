/*
 *  Created by maruyama on 11/04/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 */

#include "boid.h"
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


void Boid::set(double x, double y, double z, double vx, double vy, double vz)
{
	position.x = x;
	position.y = y;
	position.z = z;
	velocity.x = x;
	velocity.y = y;
	velocity.z = z;
}


void Boid::get(double* x, double* y, double* z, double* vx, double* vy, double* vz)
{
	*x = position.x;
	*y = position.y;
	*z = position.z;
	*vx = velocity.x;
	*vy = velocity.y;
	*vz = velocity.z;
}


void Boid::set_serialized_data(double* buffer)
{
	position.x = buffer[0];
	position.y = buffer[1];
	position.z = buffer[2];
	velocity.x = buffer[3];
	velocity.y = buffer[4];
	velocity.z = buffer[5];
}

void Boid::get_serialized_data(double* buffer)
{
	buffer[0] = position.x;
	buffer[1] = position.y;
	buffer[2] = position.z;
	buffer[3] = velocity.x;
	buffer[4] = velocity.y;
	buffer[5] = velocity.z;
}
/*
bool Boid::isInsideSeparationArea(Boid &_b) {
    return isInsideArea(_b, SIGHT_DISTANCE_SEPARATION, SIGHT_ANGLE_SEPARATION);
}

bool Boid::isInsideAlignmentArea(Boid &_b) {
    return isInsideArea(_b, SIGHT_DISTANCE_ALIGNMENT, SIGHT_ANGLE_ALIGNMENT);
}

bool Boid::isInsideCohesionArea(Boid &_b) {
    return isInsideArea(_b, SIGHT_DISTANCE_COHESION, SIGHT_ANGLE_COHESION);
}
*/
