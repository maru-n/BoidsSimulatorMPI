/*
 *  vector3D.h
 *  BoidTest
 *
 *  Created by maruyama on 11/04/26.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef parameter_H
#include "parameter.h"
#endif

#include <math.h>

class Vector3D{
public:
	Vector3D();
	Vector3D(double _x, double _y, double _z);
	//Vector3D& Vector3D::operator=(const Vector3D& v);
	Vector3D& operator=(const Vector3D& v);
	Vector3D& operator+=(const Vector3D& v);
	Vector3D& operator-=(const Vector3D& v);
	Vector3D& operator*=(double k);
	Vector3D& operator/=(double k);
	Vector3D operator+();
	Vector3D operator-();
	
	double norm();
	Vector3D normalized();

	double x;
	double y;
	double z;
};


//二項演算子
Vector3D operator+(const Vector3D& u,const Vector3D& v);

Vector3D operator-(const Vector3D& u,const Vector3D& v);

double operator*(const Vector3D& u,const Vector3D& v);

Vector3D operator*(const Vector3D& v, double k);
Vector3D operator*(double k ,const Vector3D& v);
Vector3D operator/(const Vector3D& v, double k);

//画面への出力
#include <iostream>
std::ostream& operator<<(std::ostream& stream, const Vector3D& v);
