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

using namespace std;

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

    this->INTERACTION_RANGE = (float) max(max(this->cohesion.sight_distance, this->alignment.sight_distance), this->separation.sight_distance);
    this->GRID_NUM = int(field_size / (INTERACTION_RANGE*2));
    this->grid.resize(GRID_NUM*GRID_NUM*GRID_NUM);
    this->GRID_SIZE = float(this->field_size / this->GRID_NUM);


    if (INTERACTION_RANGE * 2 >= this->GRID_SIZE) {
        std::cerr << "Invalid GRID_NUM!!!" << std::endl;
        exit(-1);
    }

    cout << "# Grid List Optimization enabled" << endl
         << "  grid num:"<< this->GRID_NUM << endl;
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
            unsigned int grid_x = int(boids[i].position.x / field_size_X);
            unsigned int grid_y = int(boids[i].position.x / field_size_X);
            unsigned int grid_z = int(boids[i].position.x / field_size_X);

            update_list_grid(boids[i]);
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

            update_list_grid(boids[i]);
        }
    } else {
    }
}



void BoidSimulation::update()
{
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

        Boid boid_this = boids[i];



        std::list<Boid*> target_boids_list = this->grid[boids[i].grid_index[0]];

        std::list<Boid*>::iterator itr = target_boids_list.begin();
        //for (int j = 0; j < target_boids_list.size(); ++j) {
        //Boid* boid_that = target_boids_list.inde
        for(; itr != target_boids_list.end(); ++itr) {
            Boid* boid_that = *itr;

        //for(auto boid_that = target_boids_list.begin(); boid_that != target_boids_list.end(); ++boid_that) {
            //Vector3D boids_j_pos_tmp =  (*boid_that)->position;
            //Vector3D boids_j_pos_tmp =  boid_that->position;
            Vector3D boids_j_pos_tmp =  (*itr)->position;


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
            if( boids[i].id != boid_that->id ){
            //if( i != j ){
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
                    //dv_ali[i] += boids[j].velocity;
                    //dv_ali[i] += target_boids_list[j].velocity;
                    dv_ali[i] += boid_that->velocity;
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

    clear_list_grid();

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

        update_list_grid(boids[i]);
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
           << "#OpenMP: ";
    if (boidsim.is_openmp_enabled()) {
        stream << "Enabled (max threads = " << boidsim.get_max_threads() << ")" << std::endl;
    } else {
        stream << "Disabled" << std::endl;
    }
    return stream;
}

void BoidSimulation::clear_list_grid() {
    for (int i = 0; i < grid.size(); ++i) {
        grid[i].clear();
    }
}

void BoidSimulation::update_list_grid(Boid &boid) {
    /*
    for (int i = 0; i < 8; ++i) {
        if (boid.grid_index[i] < 0) {
            continue;
        }
        this->grid[boid.grid_index[i]].remove(&boid);
    }*/

    // update boid grid info
    // boid.grid_index (int[8])
    int grid_index[3][2];
    grid_index[0][0] = int(boid.position.x / this->GRID_SIZE);
    grid_index[1][0] = int(boid.position.y / this->GRID_SIZE);
    grid_index[2][0] = int(boid.position.z / this->GRID_SIZE);

    float pos_in_grid[3];
    pos_in_grid[0] = boid.position.x - grid_index[0][0] * this->GRID_SIZE;
    pos_in_grid[1] = boid.position.y - grid_index[1][0] * this->GRID_SIZE;
    pos_in_grid[2] = boid.position.z - grid_index[2][0] * this->GRID_SIZE;

    for (int i = 0; i <3; ++i) {
        if (pos_in_grid[i] < INTERACTION_RANGE) {
            grid_index[i][1] = (grid_index[i][0] - 1 + GRID_NUM) % GRID_NUM;
        } else if (pos_in_grid[i] >= this->GRID_SIZE - INTERACTION_RANGE) {
            grid_index[i][1] = (grid_index[i][0] + 1) % GRID_NUM;
        } else {
            grid_index[i][1] = -1;
        }
    }

    // *** boid.grid_indes ***
    // first index is primally grid
    // other indices are overlapped grid and -1 mean no overlapping
    unsigned int n = 0;
    for (int i = 0; i < 2; ++i) {
        if (grid_index[0][i] == -1) continue;
        for (int j = 0; j < 2; ++j) {
            if (grid_index[1][j] == -1) continue;
            for (int k = 0; k < 2; ++k) {
                if (grid_index[2][k] == -1) continue;
                boid.grid_index[n] = grid_index[0][i] + grid_index[1][j] * this->GRID_NUM + grid_index[2][k] * this->GRID_NUM * this->GRID_NUM;
                n++;
            }
        }
    }

    while (n < 8) {
        boid.grid_index[n] = -1;
        n++;
    }

    for (int i = 0; i < 8; ++i) {
        if (boid.grid_index[i] < 0) {
            continue;
        }
        this->grid[boid.grid_index[i]].push_back(&boid);
    }
}


