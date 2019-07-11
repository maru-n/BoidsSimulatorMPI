//
// Created by Maruyama Norihiro on 2016/11/18.
//

#ifndef MASSIVESWARM_BOID_SIMULATION_H
#define MASSIVESWARM_BOID_SIMULATION_H

#include "boid.h"
#include "dtype.h"
#include "vector3D.h"
#include <vector>
#include <list>


class BoidSimulation {
public:
    BoidSimulation();
    virtual ~BoidSimulation();
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
               double velocity_min,
               std::string initialization,
               int rand_seed);
    virtual void init();
    virtual void update();
    //virtual int get(unsigned int id, double* x, double* y, double* z);
    virtual int get(unsigned int id, float* x, float* y, float* z);
    virtual int get_v(unsigned int id, float* vx, float* vy, float* vz);
    virtual int get(unsigned int id, float* x, float* y, float* z, float* vx, float* vy, float* vz);
    bool is_openmp_enabled() const;
    int get_max_threads() const;

    int get_force(unsigned int id, float* coh_x, float* coh_y, float* coh_z, float* sep_x, float* sep_y, float* sep_z, float* ali_x, float* ali_y, float* ali_z);

    unsigned int N;
    unsigned int time_step;
    double field_size;
    double field_size_X, field_size_Y, field_size_Z;
    interaction_parameters separation, alignment, cohesion;
    velocity_parameters velocity;
    Boid* boids;


    unsigned int GRID_NUM = 9;  //TODO:

protected:
    std::string initialization_type;
    int rand_seed;
    Vector3D *dv;
    Vector3D *dv_coh;
    Vector3D *dv_sep;
    Vector3D *dv_ali;

    float INTERACTION_RANGE;
    //unsigned int GRID_NUM = 9;  //TODO:
    float GRID_SIZE;
    std::vector<std::list<Boid*>> grid;

    void clear_list_grid();
    void update_list_grid(Boid &boid);
};

std::ostream& operator<<(std::ostream& stream, const BoidSimulation& boidsim);

#endif //MASSIVESWARM_BOID_SIMULATION_H
