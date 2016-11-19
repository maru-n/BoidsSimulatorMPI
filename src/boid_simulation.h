//
// Created by Maruyama Norihiro on 2016/11/18.
//

#ifndef MASSIVESWARM_BOID_SIMULATION_H
#define MASSIVESWARM_BOID_SIMULATION_H

#include "boid.h"
#include "vector3D.h"

typedef struct {
    double sight_distance;
    double sight_agnle;
    double force_coefficient;
} interaction_parameters;


class BoidSimulation {
public:
    BoidSimulation();
    void setup(unsigned int number_of_agents,
               double field_size,
               double separation_sight_distance,
               double separation_sight_angle,
               double separation_force_coefficient,
               double alignment_sight_distance,
               double alignment_sight_angle,
               double alignment_force_coefficient,
               double cohesion_sight_distance,
               double cohesion_sight_angle,
               double cohesion_force_coefficient,
               double velocity_max,
               double velocity_min);
    void update();
    bool is_openmp_enabled();
    int get_max_threads();

    unsigned int N;
    double field_size;
    interaction_parameters separation;
    interaction_parameters alignment;
    interaction_parameters cohesion;
    struct {
        double max;
        double min;
    } velocity;
    Boid* boids;

private:
    Vector3D *dv;
    Vector3D *dv_coh;
    Vector3D *dv_sep;
    Vector3D *dv_ali;
};


#endif //MASSIVESWARM_BOID_SIMULATION_H
