//
// Created by Maruyama Norihiro on 2016/11/18.
//

#include "boid_simulation.h"
#include <math.h>
#include <stdlib.h>

#define ENABLE_OPENMP

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#include <omp.h>
#endif

BoidSimulation::BoidSimulation() {}


void BoidSimulation::setup(unsigned int number_of_agents,
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
                           double velocity_min) {
    this->field_size = field_size;
    separation = {separation_sight_distance, separation_sight_angle, separation_force_coefficient};
    alignment = {alignment_sight_distance, alignment_sight_angle, alignment_force_coefficient};
    cohesion = {cohesion_sight_distance, cohesion_sight_angle, cohesion_force_coefficient};
    velocity = {velocity_max, velocity_min};
    N = number_of_agents;
    boids = new Boid[N];
    dv = new Vector3D[N];
    dv_ali = new Vector3D[N];
    dv_coh = new Vector3D[N];
    dv_sep = new Vector3D[N];
}

void BoidSimulation::init(std::string init_condition)
{

    if (init_condition == "test") {
        srand(12345);
        for (int i = 0; i < N; i++) {
            boids[i].position.x = field_size * 3 / 8 + drand48() * field_size / 4;
            boids[i].position.y = field_size * 3 / 8 + drand48() * field_size / 4;
            if (drand48() > 0.5) {
                boids[i].position.z = drand48() * field_size / 4 + field_size * 3 / 4;
            } else {
                boids[i].position.z = drand48() * field_size / 4;
            }

            double v = drand48() * (velocity.max - velocity.min) + velocity.min;
            double th1 = drand48() * M_PI;
            double th2 = drand48() * 2.0 * M_PI;
            boids[i].velocity.x = double(v * sin(th1) * cos(th2));
            boids[i].velocity.y = double(v * sin(th1) * sin(th2));
            boids[i].velocity.z = double(v * cos(th1));
        }
    } else if (init_condition=="random_uniform") {

    }
}



void BoidSimulation::update()
{
#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp parallel
#endif

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp for
#endif
    for(int i=0; i<N; i++){
        dv_coh[i].x = dv_coh[i].y = dv_coh[i].z =
        dv_sep[i].x = dv_sep[i].y = dv_sep[i].z =
        dv_ali[i].x = dv_ali[i].y = dv_ali[i].z = 0.0;
        int neivers_num_coh = 0;
        int neivers_num_sep = 0;
        int neivers_num_ali = 0;
        for(int j=0; j<N; j++){
            Vector3D boids_j_pos_tmp = boids[j].position;

            if ((boids_j_pos_tmp.x - boids[i].position.x) > field_size/2) {
                boids_j_pos_tmp.x -= field_size;
            } else if (boids[i].position.x - boids_j_pos_tmp.x > field_size/2) {
                boids_j_pos_tmp.x += field_size;
            }
            if ((boids_j_pos_tmp.y - boids[i].position.y) > field_size/2) {
                boids_j_pos_tmp.y -= field_size;
            } else if (boids[i].position.y - boids_j_pos_tmp.y > field_size/2) {
                boids_j_pos_tmp.y += field_size;
            }
            if ((boids_j_pos_tmp.z - boids[i].position.z) > field_size/2) {
                boids_j_pos_tmp.z -= field_size;
            } else if (boids[i].position.z - boids_j_pos_tmp.z > field_size/2) {
                boids_j_pos_tmp.z += field_size;
            }

            Boid target_boid(boids_j_pos_tmp);
            if( i != j ){
                // Cohesion
                if (boids[i].isInsideArea(target_boid, cohesion.sight_distance, cohesion.sight_agnle)){
                    neivers_num_coh ++;
                    dv_coh[i] += target_boid.position;
                }
                // Separation
                if (boids[i].isInsideArea(target_boid, separation.sight_distance, separation.sight_agnle)) {
                    neivers_num_sep ++;
                    dv_sep[i] += (boids[i].position - target_boid.position).normalized();
                }
                // Alignment
                if (boids[i].isInsideArea(target_boid, alignment.sight_distance, alignment.sight_agnle)) {
                    neivers_num_ali ++;
                    dv_ali[i] += boids[j].velocity;
                }
            }
        }
        if (neivers_num_coh != 0) {
            dv_coh[i] = dv_coh[i] / neivers_num_coh - boids[i].position;
        }
        if (neivers_num_sep != 0) {
            //dv_sep[i] = dv_sep[i] / neivers_num_sep;
        }
        if (neivers_num_ali != 0) {
            dv_ali[i] = dv_ali[i] / neivers_num_ali - boids[i].velocity;
        }
        dv[i] = cohesion.force_coefficient*dv_coh[i] + separation.force_coefficient*dv_sep[i] + alignment.force_coefficient*dv_ali[i];
    }

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp for
#endif
    for(int i=0; i<N; i++) {

        boids[i].velocity += dv[i];

        if(boids[i].velocity.norm()>0. && boids[i].velocity.norm()>velocity.max){
            boids[i].velocity = boids[i].velocity.normalized() * velocity.max;
        }else if(boids[i].velocity.norm()>0. && boids[i].velocity.norm()<velocity.min){
            boids[i].velocity = boids[i].velocity.normalized() * velocity.min;
        }

        //update boid
        boids[i].position += boids[i].velocity;

        //Boundary conditon
        if(boids[i].position.x < 0.0) {
            boids[i].position.x = field_size + boids[i].position.x;
        }
        if(boids[i].position.y < 0.0) {
            boids[i].position.y = field_size + boids[i].position.y;
        }
        if(boids[i].position.z < 0.0) {
            boids[i].position.z = field_size + boids[i].position.z;
        }
        if(boids[i].position.x > field_size) {
            boids[i].position.x = boids[i].position.x - field_size;
        }
        if(boids[i].position.y > field_size) {
            boids[i].position.y = boids[i].position.y - field_size;
        }
        if(boids[i].position.z > field_size) {
            boids[i].position.z = boids[i].position.z - field_size;
        }
    }
}

bool BoidSimulation::is_openmp_enabled()
{
#if defined(_OPENMP) && defined(ENABLE_OPENMP)
    return true;
#else
    return false;
#endif
}

int BoidSimulation::get_max_threads()
{
#if defined(_OPENMP) && defined(ENABLE_OPENMP)
    return omp_get_max_threads();
#else
    return -1;
#endif
}