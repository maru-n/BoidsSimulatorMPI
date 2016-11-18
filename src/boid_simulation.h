//
// Created by Maruyama Norihiro on 2016/11/18.
//

#ifndef MASSIVESWARM_BOID_SIMULATION_H
#define MASSIVESWARM_BOID_SIMULATION_H

#include "boid.h"
#include "vector3D.h"

class BoidSimulation {
public:
    BoidSimulation();
    void init(unsigned int number_of_agents);
    void update();
    bool is_openmp_enabled();
    int get_max_threads();

    Boid* boids;

private:
    unsigned int N;
    Vector3D *dv;
    Vector3D *dv_coh;
    Vector3D *dv_sep;
    Vector3D *dv_ali;
};


#endif //MASSIVESWARM_BOID_SIMULATION_H
