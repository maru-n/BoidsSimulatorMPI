/*
 *  Created by maruyama on 11/04/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 */

#include "boid.h"
#include <math.h>
static const double PI = 6*asin( 0.5 );

Boid::Boid(){
	Vector3D p;
	Vector3D v;
	position = p;
	velocity = v;
	live = false;
}

Boid::Boid(Vector3D _p){
	Vector3D v;
	position = _p;
	velocity = v;
	live = false;
}

Boid::Boid(Vector3D _p, Vector3D _v){
	position = _p;
	velocity = _v;
	live = false;
}
/*
void Boid::draw() {

	double bodySize = 0.01*sqrt(5.0+2*sqrt(5.0));
	double wingSize = 0.005;
	bool f=false;
    if(velocity.x==0.0 or velocity.y==0.0){
		f=true;
	}
	if(!f){
        Vector3D head = velocity.getUnity();
        head = head*bodySize;
        Vector3D wing1( wingSize * head.y / sqrt( head.x*head.x+head.y*head.y),
                       -wingSize * head.x / sqrt( head.x*head.x+head.y*head.y),
                       0.0);
        Vector3D wing2 = -wing1;
        head = position + head;
        wing1 = position + wing1;
        wing2 = position + wing2;
        glBegin(GL_TRIANGLES);
        glVertex3d(head.x, head.y, head.z);
        glVertex3d(wing1.x, wing1.y, wing1.z);
        glVertex3d(wing2.x, wing2.y, wing2.z);
        glEnd();
    }else{
        glBegin(GL_TRIANGLES);
        for (int i=0; i<10; i++) {
            double ang = 2*PI*i/10.0;
            glVertex3d(position.x+sin(ang)*bodySize, position.y+cos(ang)*bodySize, position.z);
        }
        glEnd();
    }

    //glBegin(GL_POINT);
    //glVertex3d(position.x, position.y, position.z);
    //glEnd();
}
*/
bool Boid::isInsideArea(Boid _b, double _sighe_distance, double _sight_angle) {
	if( ( (*this).position - _b.position ).getAbs() < _sighe_distance ){
        if (_sight_angle >= 1.0) {
            return true;
        }
		Vector3D sightVec = (*this).velocity;
        Vector3D relativePos = _b.position - (*this).position;
        double th = acos((sightVec * relativePos) / (sightVec.getAbs() * relativePos.getAbs()));
		if( th <= PI*_sight_angle ){
			return true;
        }else{
            return false;
        }
	}
	return false;
}

bool Boid::isInsideSeparationArea(Boid _b) {
    return isInsideArea(_b, SIGHT_DISTANCE_SEPARATION, SIGHT_ANGLE_SEPARATION);
}

bool Boid::isInsideAlignmentArea(Boid _b) {
    return isInsideArea(_b, SIGHT_DISTANCE_ALIGNMENT, SIGHT_ANGLE_ALIGNMENT);
}

bool Boid::isInsideCohesionArea(Boid _b) {
    return isInsideArea(_b, SIGHT_DISTANCE_COHESION, SIGHT_ANGLE_COHESION);
}

