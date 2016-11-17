/*
 *  parameter.h
 *  BoidTest
 *
 *  Created by maruyama on 11/04/26.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef parameter_H
#define parameter_H
#include <string>

const unsigned int fps = 30;

const double FIELD_SIZE = 1;

const double SIGHT_DISTANCE_SEPARATION = 0.1;
const double SIGHT_DISTANCE_ALIGNMENT = 0.1;
const double SIGHT_DISTANCE_COHESION = 0.1;
const double SIGHT_ANGLE_SEPARATION = 1; // 09
const double SIGHT_ANGLE_ALIGNMENT = 1; // 07
const double SIGHT_ANGLE_COHESION = 1; //08

const double COEFF_SEPARATION = 0.0001; // 0.00002
const double COEFF_COHESION = 0.1; // 0.02
const double COEFF_ALIGNMENT = 0.1; // 0.1

const double MIN_VELOCITY = 0.001;
const double MAX_VELOCITY = 0.003;


#endif
