/*
 *  vector3D.h
 *  BoidTest
 *
 *  Created by maruyama on 11/04/26.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __BOID_VECTOR_3D_H__
#define __BOID_VECTOR_3D_H__

class Vector3D{
public:
	Vector3D();
	Vector3D(float _x, float _y, float _z);
	//Vector3D& Vector3D::operator=(const Vector3D& v);
	Vector3D& operator=(const Vector3D& v);
	Vector3D& operator+=(const Vector3D& v);
	Vector3D& operator-=(const Vector3D& v);
	Vector3D& operator*=(float k);
	Vector3D& operator/=(float k);
	Vector3D operator+();
	Vector3D operator-();
	
	float norm();
	Vector3D normalized();

    float getAbs(){return norm();};  // legacy
    Vector3D getUnity(){return normalized();};  // legacy

	//double x, y, z;
	float x, y, z;
};


Vector3D operator+(const Vector3D& u,const Vector3D& v);

Vector3D operator-(const Vector3D& u,const Vector3D& v);

float operator*(const Vector3D& u,const Vector3D& v);

Vector3D operator*(const Vector3D& v, float k);
Vector3D operator*(float k ,const Vector3D& v);
Vector3D operator/(const Vector3D& v, float k);

#include <iostream>
std::ostream& operator<<(std::ostream& stream, const Vector3D& v);

#endif
