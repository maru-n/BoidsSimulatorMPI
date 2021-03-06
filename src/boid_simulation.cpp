//
// Created by Maruyama Norihiro on 2016/11/18.
//

#include "boid_simulation.h"
#include "dtype.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifdef MPI_ENABLE
#include "mpi.h"
#include <mpi-ext.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif


template <typename T> T fix_byte_order2(T value) {
    int x = 1;   // 0x00000001
    if (*(char *)&x != 0) {
        return value;  //is little endian
    }
    T ret;
    int size = sizeof(T);
    for (int i = 0; i < size; ++i) {
        ((char*)&ret)[i] = ((char*)&value)[size-i-1];
    }
    return ret;
}

BoidSimulation::BoidSimulation() {}

BoidSimulation::~BoidSimulation()
{
    delete boids;
    delete dv;
    delete dv_coh;
    delete dv_sep;
    delete dv_ali;
}

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
                           double velocity_min,
                           std::string initialization,
                           int rand_seed)
{
    field_size_X = field_size_Y = field_size_Z = field_size;
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
    this->initialization_type = initialization;
    this->rand_seed = rand_seed;
    this->time_step = 0;
}

void BoidSimulation::init()
{
    srand(rand_seed);
    if (initialization_type == "test") {
        for (int i = 0; i < N; i++) {
            boids[i].position.x = field_size_X * 3 / 8 + drand48() * field_size_X / 4;
            boids[i].position.y = field_size_Y * 3 / 8 + drand48() * field_size_Y / 4;
            if (drand48() > 0.5) {
                boids[i].position.z = drand48() * field_size_Z / 4 + field_size_Z * 3 / 4;
            } else {
                boids[i].position.z = drand48() * field_size_Z / 4;
            }

            double v = drand48() * (velocity.max - velocity.min) + velocity.min;
            double th1 = drand48() * M_PI;
            double th2 = drand48() * 2.0 * M_PI;
            boids[i].velocity.x = double(v * sin(th1) * cos(th2));
            boids[i].velocity.y = double(v * sin(th1) * sin(th2));
            boids[i].velocity.z = double(v * cos(th1));
            boids[i].id = i;
        }
    } else if (initialization_type == "random_uniform") {
        for (int i = 0; i < N; i++) {
            boids[i].position.x = field_size_X * drand48();
            boids[i].position.y = field_size_Y * drand48();
            boids[i].position.z = field_size_Z * drand48();
            double v = drand48() * (velocity.max - velocity.min) + velocity.min;
            double th1 = drand48() * M_PI;
            double th2 = drand48() * 2.0 * M_PI;
            boids[i].velocity.x = double(v * sin(th1) * cos(th2));
            boids[i].velocity.y = double(v * sin(th1) * sin(th2));
            boids[i].velocity.z = double(v * cos(th1));
            boids[i].id = i;
        }
    } else {
        std::ifstream init_data_file(initialization_type,  std::ios::binary);
        if (!init_data_file) {
            std::cerr << "Invalid initial placement type or file name: " << initialization_type << std::endl;
            exit(-1);
        }
        for (int i = 0; i < this->N; ++i) {
            boids[i].id = i;
            float tmp;
            init_data_file.read((char*)&tmp, sizeof(tmp));
            boids[i].position.x = fix_byte_order2(tmp);
            init_data_file.read((char*)&tmp, sizeof(tmp));
            boids[i].position.y = fix_byte_order2(tmp);
            init_data_file.read((char*)&tmp, sizeof(tmp));
            boids[i].position.z = fix_byte_order2(tmp);
            init_data_file.read((char*)&tmp, sizeof(tmp));
            boids[i].velocity.x = fix_byte_order2(tmp);
            init_data_file.read((char*)&tmp, sizeof(tmp));
            boids[i].velocity.y = fix_byte_order2(tmp);
            init_data_file.read((char*)&tmp, sizeof(tmp));
            boids[i].velocity.z = fix_byte_order2(tmp);
        }
        init_data_file.close();
    }
}

#include <iomanip>
#include <limits>


void BoidSimulation::update()
{
    //std::cerr << std::setprecision(std::numeric_limits<float>::max_digits10) << boids[0].position.x << "," << boids[0].velocity.x << std::endl;
#ifdef _OPENMP
#pragma omp parallel
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for(int i=0; i<N; i++){
        dv_coh[i].x = dv_coh[i].y = dv_coh[i].z =
        dv_sep[i].x = dv_sep[i].y = dv_sep[i].z =
        dv_ali[i].x = dv_ali[i].y = dv_ali[i].z = 0.0;
        int neivers_num_coh = 0;
        int neivers_num_sep = 0;
        int neivers_num_ali = 0;
        //if(i==13) {std::cerr << "#13 see ";}
        for(int j=0; j<N; j++){
            Vector3D boids_j_pos_tmp = boids[j].position;

            if ((boids_j_pos_tmp.x - boids[i].position.x) > field_size_X/2) {
                boids_j_pos_tmp.x -= field_size_X;
            } else if (boids[i].position.x - boids_j_pos_tmp.x > field_size_X/2) {
                boids_j_pos_tmp.x += field_size_X;
            }
            if ((boids_j_pos_tmp.y - boids[i].position.y) > field_size_Y/2) {
                boids_j_pos_tmp.y -= field_size_Y;
            } else if (boids[i].position.y - boids_j_pos_tmp.y > field_size_Y/2) {
                boids_j_pos_tmp.y += field_size_Y;
            }
            if ((boids_j_pos_tmp.z - boids[i].position.z) > field_size_Z/2) {
                boids_j_pos_tmp.z -= field_size_Z;
            } else if (boids[i].position.z - boids_j_pos_tmp.z > field_size_Z/2) {
                boids_j_pos_tmp.z += field_size_Z;
            }

            Boid target_boid(boids_j_pos_tmp);
            if( i != j ){
                // Cohesion
                if (boids[i].isInsideArea(target_boid, cohesion.sight_distance, cohesion.sight_agnle)){
                    //if(i==13) {std::cerr << j << ", ";}
                    neivers_num_coh ++;
                    dv_coh[i] += target_boid.position;
                }
                // Separation
                if (boids[i].isInsideArea(target_boid, separation.sight_distance, separation.sight_agnle)) {
                    //if(i==13) {std::cerr << j << ", ";}
                    neivers_num_sep ++;
                    dv_sep[i] += (boids[i].position - target_boid.position).normalized();
                }
                // Alignment
                if (boids[i].isInsideArea(target_boid, alignment.sight_distance, alignment.sight_agnle)) {
                    //if(i==13) {std::cerr << j << ", ";}
                    neivers_num_ali ++;
                    dv_ali[i] += boids[j].velocity;
                }
            }
        }
        //if(i==13) {std::cerr << std::endl;}
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

#ifdef _OPENMP
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
            boids[i].position.x = field_size_X + boids[i].position.x;

        } else if(boids[i].position.x > field_size_X) {
            boids[i].position.x = boids[i].position.x - field_size_X;
        }
        if(boids[i].position.y < 0.0) {
            boids[i].position.y = field_size_Y + boids[i].position.y;
        } else if(boids[i].position.y > field_size_Y) {
            boids[i].position.y = boids[i].position.y - field_size_Y;
        }
        if(boids[i].position.z < 0.0) {
            boids[i].position.z = field_size_Z + boids[i].position.z;
        } else if(boids[i].position.z > field_size_Z) {
            boids[i].position.z = boids[i].position.z - field_size_Z;
        }
    }
}

//int BoidSimulation::get(unsigned int id, double* x, double* y, double* z)
int BoidSimulation::get(unsigned int id, float* x, float* y, float* z)
{
    *x = boids[id].position.x;
    *y = boids[id].position.y;
    *z = boids[id].position.z;
    return 0;
}

int BoidSimulation::get_v(unsigned int id, float* vx, float* vy, float* vz)
{
    *vx = boids[id].velocity.x;
    *vy = boids[id].velocity.y;
    *vz = boids[id].velocity.z;
    return 0;
}

int BoidSimulation::get(unsigned int id, float* x, float* y, float* z, float* vx, float* vy, float* vz)
{
    this->get(id, x , y, z);
    this->get_v(id, vx, vy, vz);
    return 0;
}

int BoidSimulation::get_force(unsigned int id, float* coh_x, float* coh_y, float* coh_z, float* sep_x, float* sep_y, float* sep_z, float* ali_x, float* ali_y, float* ali_z)
{
    *coh_x = (float)cohesion.force_coefficient*dv_coh[id].x;
    *coh_y = (float)cohesion.force_coefficient*dv_coh[id].y;
    *coh_z = (float)cohesion.force_coefficient*dv_coh[id].z;

    *sep_x = (float)separation.force_coefficient*dv_sep[id].x;
    *sep_y = (float)separation.force_coefficient*dv_sep[id].y;
    *sep_z = (float)separation.force_coefficient*dv_sep[id].z;

    *ali_x = (float)alignment.force_coefficient*dv_ali[id].x;
    *ali_y = (float)alignment.force_coefficient*dv_ali[id].y;
    *ali_z = (float)alignment.force_coefficient*dv_ali[id].z;

    return 0;
}


bool BoidSimulation::is_openmp_enabled() const
{
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

int BoidSimulation::get_max_threads() const
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return -1;
#endif
}

std::ostream& operator<<(std::ostream& stream, const BoidSimulation& boidsim)
{
    stream << "N: " << boidsim.N << std::endl
           << "field size: " << boidsim.field_size_X << "," << boidsim.field_size_Y << "," << boidsim.field_size_Z << std::endl
           << "#Separation" << std::endl
           << "force: " << boidsim.separation.force_coefficient << std::endl
           << "area distance: " << boidsim.separation.sight_distance << std::endl
           << "area angle: " << boidsim.separation.sight_agnle << std::endl
           << "#Alignment" << std::endl
           << "force: " << boidsim.alignment.force_coefficient << std::endl
           << "area distance: " <<  boidsim.alignment.sight_distance << std::endl
           << "area angle: " << boidsim.alignment.sight_agnle << std::endl
           << "#Cohesion" << std::endl
           << "force: " << boidsim.cohesion.force_coefficient << std::endl
           << "area distance: " << boidsim.cohesion.sight_distance << std::endl
           << "area angle: " << boidsim.cohesion.sight_agnle << std::endl
           << "#Velocity" << std::endl
           << "min: " << boidsim.velocity.min << std::endl
           << "max: " << boidsim.velocity.max << std::endl
           << "#Grid List Optimization: " << boidsim.GRID_NUM << std::endl
           << "#OpenMP: ";
    if (boidsim.is_openmp_enabled()) {
        stream << "Enabled (max threads = " << boidsim.get_max_threads() << ")" << std::endl;
    } else {
        stream << "Disabled" << std::endl;
    }
    return stream;
}
