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

//const std::string fname = "/Users/maruyama/Desktop/data_local.ptcl";
//const std::string fname = "/Users/maruyama/Desktop/data.ptcl";
//const std::string fname = "~/data.ptcl";
//const std::string fname = "/data/hp160264/k03378/data.ptcl";
//const std::string fname = "/data/hp160264/k03378/data.ptcl";
//const std::string fname = "/Users/maruyama/Desktop/data_local.ptcl";
//const std::string fname = "data.ptcl";
//std::string fname = "/Users/maruyama/Desktop/testdata.ptcl";
//std::string fname;

//const int N = 256;
//const int N =256;
//unsigned int N;
//const unsigned int N = 1000;
//const unsigned int N = 256;
//const unsigned int N = 10000;
//unsigned int N;
//unsigned int T;
//const unsigned int T = 1800;
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

//const double preyF = 0.00000;

//const int toLiveMin = 10; //if( toLiveMin <= seeing boids <= toLiveMax )
//const int toLiveMax = 25; //the boid tweet

const double MIN_VELOCITY = 0.001;
const double MAX_VELOCITY = 0.003;

//const bool LOGGING = true;
//const int CHANGE_PREY = 3000; //step
//const int outputStep = 1;
//const int exitStep = 5000; //step, when only the next param is true, it's work.

//const bool stopPeriod = false;
//const bool PREY = false;
//const bool realTime = false; //don't use this

#endif
