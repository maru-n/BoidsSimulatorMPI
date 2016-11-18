//
// Created by Maruyama Norihiro on 2016/11/18.
//

#include "boid_simulation.h"
#include <math.h>
#include <stdlib.h>
#include "parameter.h"

#define ENABLE_OPENMP

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#include <omp.h>
#endif

BoidSimulation::BoidSimulation()
{
}

void BoidSimulation::init(unsigned int number_of_agents)
{
    N = number_of_agents;
    boids = new Boid[N];
    dv = new Vector3D[N];
    dv_ali = new Vector3D[N];
    dv_coh = new Vector3D[N];
    dv_sep = new Vector3D[N];

    srand(12345);
    for(int i=0; i<N; i++){
        boids[i].position.x =FIELD_SIZE*3/8 + drand48()*FIELD_SIZE/4;
        boids[i].position.y =FIELD_SIZE*3/8 + drand48()*FIELD_SIZE/4;
        if (drand48() > 0.5) {
            boids[i].position.z = drand48()*FIELD_SIZE/4 + FIELD_SIZE*3/4;
        } else {
            boids[i].position.z = drand48()*FIELD_SIZE/4;
        }

        double v = drand48() * (MAX_VELOCITY - MIN_VELOCITY) + MIN_VELOCITY;
        double th1 = drand48() * M_PI;
        double th2 = drand48() * 2.0 * M_PI;
        boids[i].velocity.x = double(v * sin(th1) * cos(th2));
        boids[i].velocity.y = double(v * sin(th1) * sin(th2));
        boids[i].velocity.z = double(v * cos(th1));
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

            if ((boids_j_pos_tmp.x - boids[i].position.x) > FIELD_SIZE/2) {
                boids_j_pos_tmp.x -= FIELD_SIZE;
            } else if (boids[i].position.x - boids_j_pos_tmp.x > FIELD_SIZE/2) {
                boids_j_pos_tmp.x += FIELD_SIZE;
            }
            if ((boids_j_pos_tmp.y - boids[i].position.y) > FIELD_SIZE/2) {
                boids_j_pos_tmp.y -= FIELD_SIZE;
            } else if (boids[i].position.y - boids_j_pos_tmp.y > FIELD_SIZE/2) {
                boids_j_pos_tmp.y += FIELD_SIZE;
            }
            if ((boids_j_pos_tmp.z - boids[i].position.z) > FIELD_SIZE/2) {
                boids_j_pos_tmp.z -= FIELD_SIZE;
            } else if (boids[i].position.z - boids_j_pos_tmp.z > FIELD_SIZE/2) {
                boids_j_pos_tmp.z += FIELD_SIZE;
            }

            Boid target_boid(boids_j_pos_tmp);
            if( i != j ){
                // Cohesion
                if (boids[i].isInsideCohesionArea(target_boid)){
                    neivers_num_coh ++;
                    dv_coh[i] += target_boid.position;
                }
                // Separation
                if (boids[i].isInsideSeparationArea(target_boid)) {
                    neivers_num_sep ++;
                    dv_sep[i] += (boids[i].position - target_boid.position).normalized();
                }
                // Alignment
                if (boids[i].isInsideAlignmentArea(target_boid)) {
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
        dv[i] = COEFF_COHESION*dv_coh[i] + COEFF_SEPARATION*dv_sep[i] +COEFF_ALIGNMENT*dv_ali[i];
    }

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp for
#endif
    for(int i=0; i<N; i++) {

        boids[i].velocity += dv[i];

        if(boids[i].velocity.norm()>0. && boids[i].velocity.norm()>MAX_VELOCITY){
            boids[i].velocity = boids[i].velocity.normalized() * MAX_VELOCITY;
        }else if(boids[i].velocity.norm()>0. && boids[i].velocity.norm()<MIN_VELOCITY){
            boids[i].velocity = boids[i].velocity.normalized() * MIN_VELOCITY;
        }

        //update boid
        boids[i].position += boids[i].velocity;

        //Boundary conditon
        if(boids[i].position.x < 0.0) {
            boids[i].position.x = FIELD_SIZE + boids[i].position.x;
        }
        if(boids[i].position.y < 0.0) {
            boids[i].position.y = FIELD_SIZE + boids[i].position.y;
        }
        if(boids[i].position.z < 0.0) {
            boids[i].position.z = FIELD_SIZE + boids[i].position.z;
        }
        if(boids[i].position.x > FIELD_SIZE) {
            boids[i].position.x = boids[i].position.x - FIELD_SIZE;
        }
        if(boids[i].position.y > FIELD_SIZE) {
            boids[i].position.y = boids[i].position.y - FIELD_SIZE;
        }
        if(boids[i].position.z > FIELD_SIZE) {
            boids[i].position.z = boids[i].position.z - FIELD_SIZE;
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