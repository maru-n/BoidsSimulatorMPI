        /*
 *  vector3D.cpp
 *  BoidTest
 *
 *  Created by maruyama on 11/04/26.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "vector3D.h"
#include <math.h>

Vector3D::Vector3D():x(0),y(0),z(0)
{
}

Vector3D::Vector3D(double _x, double _y, double _z):x(_x),y(_y),z(_z)
{
}
//代入演算子の定義
Vector3D& Vector3D::operator=(const Vector3D& v){
	this->x=v.x;
	this->y=v.y;
	this->z=v.z;
	return *this;
}
// +=の定義
Vector3D& Vector3D::operator+=(const Vector3D& v){
	this->x += v.x;
	this->y += v.y;
	this->z += v.z;
	return *this;		
}
// -=の定義
Vector3D& Vector3D::operator-=(const Vector3D& v){
	this->x -= v.x;
	this->y -= v.y;
	this->z -= v.z;
	return *this;		
}
// *=の定義
Vector3D& Vector3D::operator*=(double k){
	this->x *= k;
	this->y *= k;
	this->z *= k;
	return *this;	
}
// /=の定義
Vector3D& Vector3D::operator/=(double k){
	this->x /= k;
	this->y /= k;
	this->z /= k;
	return *this;	
}
//+の定義:	+v
Vector3D Vector3D::operator+(){			
	return *this;
}
//-の定義:	-v
Vector3D Vector3D::operator-(){			
	return Vector3D(-x,-y,-z);
}

double Vector3D::norm(){
	double d = x*x + y*y + z*z;
	return sqrt(d);
}

Vector3D Vector3D::normalized(){
	Vector3D v(x,y,z);
	//Vector3D v(1,1,1);
	v /= v.norm();
	return v;
}


//二項演算子
Vector3D operator+(const Vector3D& u,const Vector3D& v){	//vector+vector
	Vector3D w;
	w.x = u.x + v.x;
	w.y = u.y + v.y;
	w.z = u.z + v.z;
	return w;
}

Vector3D operator-(const Vector3D& u,const Vector3D& v){	//vector-vector
	Vector3D w;
	w.x = u.x - v.x;
	w.y = u.y - v.y;
	w.z = u.z - v.z;
	return w;
}

double operator*(const Vector3D& u,const Vector3D& v){	//内積 vector*vector
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

Vector3D operator*(const Vector3D& v, double k){	//vector*scalar
	Vector3D w;
	w.x = v.x * k;
	w.y = v.y * k;
	w.z = v.z * k;
	return w;
}
Vector3D operator*(double k ,const Vector3D& v){	//scalar*vector
	Vector3D w;
	w.x = v.x * k;
	w.y = v.y * k;
	w.z = v.z * k;
	return w;
}
Vector3D operator/(const Vector3D& v, double k){	//vector/scalar
//Vector3D operator/(double k, const Vector3D& v){	//vector/scalar
	Vector3D w;
	w.x = v.x / k;
	w.y = v.y / k;
	w.z = v.z / k;
	return w;
}

//画面への出力
#include <iostream>
using namespace std;
ostream& operator<<(ostream& stream, const Vector3D& v){
	return stream <<'('<<v.x<<","<<v.y<<","<<v.z<<')';
}

