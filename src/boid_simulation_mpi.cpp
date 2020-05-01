//
// Created by Maruyama Norihiro on 2016/11/22.
//

#include "boid.h"
#include "boid_simulation_mpi.h"
#include "mpi.h"
#include <mpi-ext.h>
#include <cfloat>
#include <cmath>
#include <algorithm>

using namespace std;

int topology_get_dimension(int *size)
{
#ifdef __FCC_VERSION
    return FJMPI_Topology_get_dimension(size);
#else
    *size = 3;
    return 0;
#endif
}

int topology_get_shape(int *x, int *y, int *z)
{
#ifdef __FCC_VERSION
    return FJMPI_Topology_get_shape(x, y, z);
#else
    //TODO: this assume that node size is fixed 8. (mpiexec -n 8)
    *x = 2;
    *y = 2;
    *z = 2;
    return 0;
#endif
}

int topology_rank2xyz(int rank, int *x, int *y, int *z)
{
#ifdef __FCC_VERSION
    return FJMPI_Topology_rank2xyz(rank, x, y, z);
#else
    int size_x, size_y, size_z;
    topology_get_shape(&size_x, &size_y, &size_z);
    *x = rank % size_x;
    *y = rank / size_x % size_y;
    *z = rank / (size_x * size_y);
    return 0;
#endif
}

int topology_xyz2rank(int x, int y, int z, int *rank)
{
#ifdef __FCC_VERSION
    return FJMPI_Topology_xyz2rank(x, y, z, rank);
#else
    int size_x, size_y, size_z;
    topology_get_shape(&size_x, &size_y, &size_z);
    *rank = size_x*size_y*z + size_x*y + x;
    return 0;
#endif
}
int step = 0;

BoidSimulationMultiNode::BoidSimulationMultiNode(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    topology_get_dimension(&mpi_dim);
    topology_get_shape(&mpi_topology_x, &mpi_topology_y, &mpi_topology_z);
    topology_rank2xyz(mpi_rank, &mpi_position_x, &mpi_position_y, &mpi_position_z);
    is_master = (mpi_rank == 0);
    if(mpi_dim != 3) {
        std::cerr << "Invalid node topology." << std::endl;
        exit(-1);
    }

    neighborhood_num = 0;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                int n_x = (mpi_position_x+mpi_topology_x+i-1)%mpi_topology_x;
                int n_y = (mpi_position_y+mpi_topology_y+j-1)%mpi_topology_y;
                int n_z = (mpi_position_z+mpi_topology_z+k-1)%mpi_topology_z;
                int r;
                topology_xyz2rank(n_x, n_y, n_z, &r);
                bool new_rank = true;
                for (int l = 0; l < neighborhood_num; ++l) {
                    if (neighborhood_rank[l] == r) {
                        new_rank = false;
                        break;
                    }
                }
                if (new_rank) {
                    neighborhood_rank[neighborhood_num++] = r;
                }
            }
        }
    }
}

BoidSimulationMultiNode::~BoidSimulationMultiNode()
{
    MPI_Finalize();
    delete[] data_buffer;
    delete[] data_id_buffer;
    delete[] data_id_buffer_swap;
    delete[] data_buffer_swap;
    delete[] margin_data_buffer;
    delete[] margin_data_buffer_swap;
    delete[] margin_data_id_buffer;
    delete[] margin_data_id_buffer_swap;
    delete[] send_data_buffer;
    delete[] send_data_id_buffer;
    delete[] recv_data_buffer;
    delete[] recv_data_id_buffer;
    delete[] data_num_buffer;
    delete[] data_id_x;
}

void BoidSimulationMultiNode::init()
{
    //data_buffer = new double[N*6];
    //data_buffer_swap = new double[N*6];
    data_buffer = new float[N*6];
    data_buffer_swap = new float[N*6];
    data_id_buffer = new unsigned[N];
    data_id_buffer_swap = new unsigned[N];

    data_id_x = new id_x[N];

    //margin_data_buffer = new double[N*6];
    //margin_data_buffer_swap = new double[N*6];
    margin_data_buffer = new float[N*6];
    margin_data_buffer_swap = new float[N*6];
    margin_data_id_buffer = new unsigned[N];
    margin_data_id_buffer_swap = new unsigned[N];

    data_num_buffer = new unsigned int[mpi_size];

    //send_data_buffer = new double[N*6];
    send_data_buffer = new float[N*6];
    send_data_id_buffer = new unsigned[N];
    //recv_data_buffer = new double[N*6];
    recv_data_buffer = new float[N*6];
    recv_data_id_buffer = new unsigned[N];

    field_size_local_x = field_size / mpi_topology_x;
    field_size_local_y = field_size / mpi_topology_y;
    field_size_local_z = field_size / mpi_topology_z;

    space_x_lower = mpi_position_x * field_size_local_x;
    space_x_upper = space_x_lower + field_size_local_x;
    space_y_lower = mpi_position_y * field_size_local_y;
    space_y_upper = space_y_lower + field_size_local_y;
    space_z_lower = mpi_position_z * field_size_local_z;
    space_z_upper = space_z_lower + field_size_local_z;

    margin_width = std::max(std::max(separation.sight_distance, alignment.sight_distance), cohesion.sight_distance);

    this->buffer_2_margin_buffer_idx = new unsigned int[N];
    field_size_local_x_with_margin = field_size_local_x + margin_width*2;
    field_size_local_y_with_margin = field_size_local_y + margin_width*2;
    field_size_local_z_with_margin = field_size_local_z + margin_width*2;

    margin_x_lower = space_x_lower - margin_width;
    if (margin_x_lower < 0.0) margin_x_lower += field_size;

    margin_y_lower = space_y_lower - margin_width;
    if (margin_y_lower < 0.0) margin_y_lower += field_size;

    margin_z_lower = space_z_lower - margin_width;
    if (margin_z_lower < 0.0) margin_z_lower += field_size;

    margin_x_upper = space_x_upper + margin_width;
    if (margin_x_upper >= field_size) margin_x_upper -= field_size;

    margin_y_upper = space_y_upper + margin_width;
    if (margin_y_upper >= field_size) margin_y_upper -= field_size;

    margin_z_upper = space_z_upper + margin_width;
    if (margin_z_upper >= field_size) margin_z_upper -= field_size;

    if (is_master) {
        BoidSimulation::init();
        for (int i = 0; i < N; ++i) {
            boids[i].get_serialized_data(&data_buffer_swap[i*6]);
        }
    }
    //MPI_Bcast(data_buffer, N*6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(data_buffer_swap, N*6, MPI_FLOAT, 0, MPI_COMM_WORLD);
    this->init_list_grid();

    unsigned id;
    //double x, y, z, vx, vy, vz;
    float x, y, z, vx, vy, vz;
    data_buffer_count = 0;
    margin_data_buffer_count = 0;
    for (int i = 0; i < N; ++i) {
        id = (unsigned)i;
        x = data_buffer_swap[i*6+0];
        y = data_buffer_swap[i*6+1];
        z = data_buffer_swap[i*6+2];
        vx = data_buffer_swap[i*6+3];
        vy = data_buffer_swap[i*6+4];
        vz = data_buffer_swap[i*6+5];

        if (x == field_size_X) x = 0.0f;
        if (y == field_size_Y) y = 0.0f;
        if (z == field_size_Z) z = 0.0f;

        /*
        if (x >= space_x_lower && (mpi_position_x+1==mpi_topology_x ? x <= space_x_upper : x < space_x_upper) &&
            y >= space_y_lower && (mpi_position_y+1==mpi_topology_y ? y <= space_y_upper : y < space_y_upper) &&
            z >= space_z_lower && (mpi_position_z+1==mpi_topology_z ? z <= space_z_upper : z < space_z_upper) ) {
            */
        if (x >= space_x_lower &&  x < space_x_upper &&
            y >= space_y_lower &&  y < space_y_upper &&
            z >= space_z_lower &&  z < space_z_upper ) {
            buffer_2_margin_buffer_idx[data_buffer_count / 6] = margin_data_buffer_count / 6;

            data_id_buffer[data_buffer_count/6] = id;
            data_buffer[data_buffer_count+0] = x;
            data_buffer[data_buffer_count+1] = y;
            data_buffer[data_buffer_count+2] = z;
            data_buffer[data_buffer_count+3] = vx;
            data_buffer[data_buffer_count+4] = vy;
            data_buffer[data_buffer_count+5] = vz;
            data_buffer_count += 6;
        }

        if ( (margin_x_lower < margin_x_upper ? (x>=margin_x_lower&&x<=margin_x_upper) : (x>=margin_x_lower||x<=margin_x_upper)) &&
             (margin_y_lower < margin_y_upper ? (y>=margin_y_lower&&y<=margin_y_upper) : (y>=margin_y_lower||y<=margin_y_upper)) &&
             (margin_z_lower < margin_z_upper ? (z>=margin_z_lower&&z<=margin_z_upper) : (z>=margin_z_lower||z<=margin_z_upper))) {
            margin_data_id_buffer[margin_data_buffer_count/6] = id;
            margin_data_buffer[margin_data_buffer_count+0] = x;
            margin_data_buffer[margin_data_buffer_count+1] = y;
            margin_data_buffer[margin_data_buffer_count+2] = z;
            margin_data_buffer[margin_data_buffer_count+3] = vx;
            margin_data_buffer[margin_data_buffer_count+4] = vy;
            margin_data_buffer[margin_data_buffer_count+5] = vz;
            margin_data_buffer_count +=6;
        }
    }

    //this->gather_data();
}


void BoidSimulationMultiNode::update()
{
    /*
     * Update
     */
    unsigned id;
    //double x, y, z, vx, vy, vz;
    float x, y, z, vx, vy, vz;


#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp parallel
#endif

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp for
//#pragma omp for schedule(guided)
#endif
    this->update_list_grid();

    for(int i=0; i<data_buffer_count/6; i++) {
        dv_coh[i].x = dv_coh[i].y = dv_coh[i].z =
        dv_sep[i].x = dv_sep[i].y = dv_sep[i].z =
        dv_ali[i].x = dv_ali[i].y = dv_ali[i].z = 0.0;
        int neivers_num_coh = 0;
        int neivers_num_sep = 0;
        int neivers_num_ali = 0;

        this_boid.set_serialized_data(&data_buffer[i*6]);
        unsigned int idx = buffer_2_margin_buffer_idx[i];

        //std::list<unsigned int> target_boids_idx_list = this->grid[grid_index[idx][0]];
        std::list<unsigned int> target_boids_idx_list = this->grid[grid_index[idx][0]];

        std::list<unsigned int>::iterator itr;
        for(itr = target_boids_idx_list.begin(); itr != target_boids_idx_list.end(); ++itr) {
            unsigned int j = *itr;
        //for(int j=0; j<margin_data_buffer_count/6; j++) {

            that_boid.set_serialized_data(&margin_data_buffer[j*6]);

            if ((that_boid.position.x - this_boid.position.x) > field_size/2) {
                that_boid.position.x -= field_size;
            } else if (this_boid.position.x - that_boid.position.x > field_size/2) {
                that_boid.position.x += field_size;
            }
            if ((that_boid.position.y - this_boid.position.y) > field_size/2) {
                that_boid.position.y -= field_size;
            } else if (this_boid.position.y - that_boid.position.y > field_size/2) {
                that_boid.position.y += field_size;
            }
            if ((that_boid.position.z - this_boid.position.z) > field_size/2) {
                that_boid.position.z -= field_size;
            } else if (this_boid.position.z - that_boid.position.z > field_size/2) {
                that_boid.position.z += field_size;
            }

            if( data_id_buffer[i] != margin_data_id_buffer[j] ){
                // Cohesion
                if (this_boid.isInsideArea(that_boid, cohesion.sight_distance, cohesion.sight_agnle)){
                    //if(data_id_buffer[i]==13) {std::cerr << margin_data_id_buffer[j] << ", ";}
                    neivers_num_coh ++;
                    dv_coh[i] += that_boid.position;
                }
                // Separation
                if (this_boid.isInsideArea(that_boid, separation.sight_distance, separation.sight_agnle)) {
                    //if(data_id_buffer[i]==13) {std::cerr << margin_data_id_buffer[j] << ", ";}
                    neivers_num_sep ++;
                    dv_sep[i] += (this_boid.position - that_boid.position).normalized();
                }
                // Alignment
                if (this_boid.isInsideArea(that_boid, alignment.sight_distance, alignment.sight_agnle)) {
                    //if(data_id_buffer[i]==13) {std::cerr << margin_data_id_buffer[j] << ", ";}
                    neivers_num_ali ++;
                    dv_ali[i] += that_boid.velocity;
                }
            }
        }
        //if(data_id_buffer[i]==13) {std::cerr << std::endl;}
        if (neivers_num_coh != 0) {
            dv_coh[i] = dv_coh[i] / neivers_num_coh - this_boid.position;
        }
        if (neivers_num_sep != 0) {
            //dv_sep[i] = dv_sep[i] / neivers_num_sep;
        }
        if (neivers_num_ali != 0) {
            dv_ali[i] = dv_ali[i] / neivers_num_ali - this_boid.velocity;
        }
        dv[i] = cohesion.force_coefficient*dv_coh[i] + separation.force_coefficient*dv_sep[i] + alignment.force_coefficient*dv_ali[i];

    }

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp for
//#pragma omp for schedule(guided)
#endif
    for(int i=0; i<data_buffer_count/6; i++) {
        boids[i].set_serialized_data(&data_buffer[i * 6]);

        //update boid
        boids[i].velocity += dv[i];
        if (boids[i].velocity.norm() > 0. && boids[i].velocity.norm() > velocity.max) {
            boids[i].velocity = boids[i].velocity.normalized() * velocity.max;
        } else if (boids[i].velocity.norm() > 0. && boids[i].velocity.norm() < velocity.min) {
            boids[i].velocity = boids[i].velocity.normalized() * velocity.min;
        }
        boids[i].position += boids[i].velocity;

        //Boundary conditon
        if (boids[i].position.x < 0.0) {
            boids[i].position.x = (float) (field_size + boids[i].position.x);
            // Take care of rounding error
            if (std::abs(boids[i].position.x - field_size) <= FLT_EPSILON ) {
                boids[i].position.x = 0.0f;
            }
        } else if (boids[i].position.x >= field_size) {
            boids[i].position.x = (float) (boids[i].position.x - field_size);
        }
        if (boids[i].position.y < 0.0) {
            boids[i].position.y = (float) (field_size + boids[i].position.y);
            // Take care of rounding error
            if (std::abs(boids[i].position.y - field_size) <= FLT_EPSILON ) {
                boids[i].position.y = 0.0f;
            }
        } else if (boids[i].position.y >= field_size) {
            boids[i].position.y = (float) (boids[i].position.y - field_size);
        }
        if (boids[i].position.z < 0.0) {
            boids[i].position.z = (float) (field_size + boids[i].position.z);
            // Take care of rounding error
            if (std::abs(boids[i].position.z - field_size) <= FLT_EPSILON ) {
                boids[i].position.z = 0.0f;
            }
        } else if (boids[i].position.z >= field_size) {
            boids[i].position.z = (float) (boids[i].position.z - field_size);
        }

        boids[i].get_serialized_data(&data_buffer[i * 6]);

        //DEBUG CODE
        /*
        unsigned bug_id = 27558;
        if (data_id_buffer[i] == bug_id) {
            std::cerr << "(" << mpi_position_x << "," << mpi_position_y << "," << mpi_position_z << ")"
                      << data_buffer[i*6+0] << ","  << data_buffer[i*6+1] << ","  << data_buffer[i*6+2] << std::endl;
        }
         */
    }


    /*
     * Comunicate with other nodes
     */

    unsigned int send_data_buffer_count = 0;
    unsigned int recv_data_buffer_count = 0;
    unsigned int data_buffer_count_new = 0;
    unsigned int margin_data_buffer_count_new = 0;

    for(int i=0; i<data_buffer_count/6; i++) {
        id = data_id_buffer[i];
        x = data_buffer[i*6+0];
        y = data_buffer[i*6+1];
        z = data_buffer[i*6+2];
        vx = data_buffer[i*6+3];
        vy = data_buffer[i*6+4];
        vz = data_buffer[i*6+5];

        if (x <= (space_x_lower+margin_width) || x >= (space_x_upper-margin_width) ||
            y <= (space_y_lower+margin_width) || y >= (space_y_upper-margin_width) ||
            z <= (space_z_lower+margin_width) || z >= (space_z_upper-margin_width)) {
            send_data_buffer[send_data_buffer_count+0] = x;
            send_data_buffer[send_data_buffer_count+1] = y;
            send_data_buffer[send_data_buffer_count+2] = z;
            send_data_buffer[send_data_buffer_count+3] = vx;
            send_data_buffer[send_data_buffer_count+4] = vy;
            send_data_buffer[send_data_buffer_count+5] = vz;
            send_data_id_buffer[send_data_buffer_count/6] = id;
            send_data_buffer_count+=6;
        } else {
            buffer_2_margin_buffer_idx[data_buffer_count_new/6] = margin_data_buffer_count_new/6;
            data_id_buffer_swap[data_buffer_count_new/6] = id;
            data_buffer_swap[data_buffer_count_new+0] = x;
            data_buffer_swap[data_buffer_count_new+1] = y;
            data_buffer_swap[data_buffer_count_new+2] = z;
            data_buffer_swap[data_buffer_count_new+3] = vx;
            data_buffer_swap[data_buffer_count_new+4] = vy;
            data_buffer_swap[data_buffer_count_new+5] = vz;
            data_buffer_count_new += 6;
            margin_data_id_buffer_swap[margin_data_buffer_count_new/6] = id;
            margin_data_buffer_swap[margin_data_buffer_count_new+0] = x;
            margin_data_buffer_swap[margin_data_buffer_count_new+1] = y;
            margin_data_buffer_swap[margin_data_buffer_count_new+2] = z;
            margin_data_buffer_swap[margin_data_buffer_count_new+3] = vx;
            margin_data_buffer_swap[margin_data_buffer_count_new+4] = vy;
            margin_data_buffer_swap[margin_data_buffer_count_new+5] = vz;
            margin_data_buffer_count_new +=6;
        }
    }
    data_buffer_count = data_buffer_count_new;
    margin_data_buffer_count = margin_data_buffer_count_new;

    //double* tmp;
    float* tmp;
    tmp = data_buffer;
    data_buffer = data_buffer_swap;
    data_buffer_swap = tmp;

    tmp = margin_data_buffer;
    margin_data_buffer = margin_data_buffer_swap;
    margin_data_buffer_swap = tmp;

    unsigned* tmp2;
    tmp2 = data_id_buffer;
    data_id_buffer = data_id_buffer_swap;
    data_id_buffer_swap = tmp2;

    tmp2 = margin_data_id_buffer;
    margin_data_id_buffer = margin_data_id_buffer_swap;
    margin_data_id_buffer_swap = tmp2;

    MPI_Request request1[27*2];
    MPI_Request request2[27*4];
    MPI_Status status1[27*2];
    MPI_Status status2[27*4];
    int request_n=0;

    unsigned int data_num[27];

    for (int i = 0; i < neighborhood_num; ++i) {
        MPI_Isend(&send_data_buffer_count, 1, MPI_INT, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request1[request_n++]);
        MPI_Irecv(&data_num[i], 1, MPI_INT, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request1[request_n++]);
    }
    MPI_Waitall(request_n, request1, status1);

    request_n = 0;
    for (int i = 0; i < neighborhood_num; ++i) {
        //MPI_Isend(send_data_buffer, send_data_buffer_count, MPI_DOUBLE, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        //MPI_Irecv(&recv_data_buffer[recv_data_buffer_count], data_num[i], MPI_DOUBLE, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Isend(send_data_buffer, send_data_buffer_count, MPI_FLOAT, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Irecv(&recv_data_buffer[recv_data_buffer_count], data_num[i], MPI_FLOAT, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Isend(send_data_id_buffer, send_data_buffer_count/6, MPI_UNSIGNED, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Irecv(&recv_data_id_buffer[recv_data_buffer_count/6], data_num[i]/6, MPI_UNSIGNED, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        recv_data_buffer_count += data_num[i];
    }

    MPI_Waitall(request_n, request2, status2);

    for (int i = 0; i < recv_data_buffer_count/6; i++) {
        id = recv_data_id_buffer[i];
        x = recv_data_buffer[i*6+0];
        y = recv_data_buffer[i*6+1];
        z = recv_data_buffer[i*6+2];
        vx = recv_data_buffer[i*6+3];
        vy = recv_data_buffer[i*6+4];
        vz = recv_data_buffer[i*6+5];
        /*
        if (x >= space_x_lower && (mpi_position_x+1==mpi_topology_x ? x <= space_x_upper : x < space_x_upper) &&
            y >= space_y_lower && (mpi_position_y+1==mpi_topology_y ? y <= space_y_upper : y < space_y_upper) &&
            z >= space_z_lower && (mpi_position_z+1==mpi_topology_z ? z <= space_z_upper : z < space_z_upper) ) {
            */
        if (x >= space_x_lower &&  x < space_x_upper &&
            y >= space_y_lower &&  y < space_y_upper &&
            z >= space_z_lower &&  z < space_z_upper ) {

            buffer_2_margin_buffer_idx[data_buffer_count/6] = margin_data_buffer_count/6;

            data_id_buffer[data_buffer_count/6] = id;
            data_buffer[data_buffer_count+0] = x;
            data_buffer[data_buffer_count+1] = y;
            data_buffer[data_buffer_count+2] = z;
            data_buffer[data_buffer_count+3] = vx;
            data_buffer[data_buffer_count+4] = vy;
            data_buffer[data_buffer_count+5] = vz;
            data_buffer_count += 6;
        }
        if ( (margin_x_lower < margin_x_upper ? (x>=margin_x_lower&&x<=margin_x_upper) : (x>=margin_x_lower||x<=margin_x_upper)) &&
             (margin_y_lower < margin_y_upper ? (y>=margin_y_lower&&y<=margin_y_upper) : (y>=margin_y_lower||y<=margin_y_upper)) &&
             (margin_z_lower < margin_z_upper ? (z>=margin_z_lower&&z<=margin_z_upper) : (z>=margin_z_lower||z<=margin_z_upper))) {
            margin_data_id_buffer[margin_data_buffer_count/6] = id;
            margin_data_buffer[margin_data_buffer_count+0] = x;
            margin_data_buffer[margin_data_buffer_count+1] = y;
            margin_data_buffer[margin_data_buffer_count+2] = z;
            margin_data_buffer[margin_data_buffer_count+3] = vx;
            margin_data_buffer[margin_data_buffer_count+4] = vy;
            margin_data_buffer[margin_data_buffer_count+5] = vz;
            margin_data_buffer_count +=6;
        }

    }
}

//int BoidSimulationMultiNode::get(unsigned int id, double* x, double* y, double* z)
int BoidSimulationMultiNode::get(unsigned int id, float* x, float* y, float* z)
{
    //assert(id == data_id_x[id].id);
    *x = data_id_x[id].x[0];
    *y = data_id_x[id].x[1];
    *z = data_id_x[id].x[2];
    return 0;
}

void BoidSimulationMultiNode::gather_data()
{
    MPI_Gather(&data_buffer_count, 1, MPI_INT, data_num_buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < data_buffer_count/6; ++i) {
        data_id_x[i].id = data_id_buffer[i];
        data_id_x[i].x[0] = data_buffer[i*6+0];
        data_id_x[i].x[1] = data_buffer[i*6+1];
        data_id_x[i].x[2] = data_buffer[i*6+2];
    }


    if (is_master) {
        MPI_Request *request = new MPI_Request[mpi_size-1];
        MPI_Status *status = new MPI_Status[mpi_size-1];
        int n = data_num_buffer[0] / 6;
        for (int i = 1; i < mpi_size; ++i) {
            MPI_Irecv(&data_id_x[n], sizeof(id_x) * data_num_buffer[i]/6, MPI_BYTE, i, i, MPI_COMM_WORLD, &request[i-1]);
            n += data_num_buffer[i] / 6;
        }
        MPI_Waitall(mpi_size-1, request, status);
        delete[] request;
        delete[] status;

        std::sort(data_id_x, data_id_x + N, [](const  id_x& x, const id_x& y) { return x.id < y.id;});
    } else {
        MPI_Request request;
        MPI_Status status;
        MPI_Isend(data_id_x, sizeof(id_x) * data_buffer_count/6, MPI_BYTE, 0, mpi_rank, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
    }
}

void BoidSimulationMultiNode::init_list_grid() {
    //this->grid_index = new int[N];
    this->grid_index = new int*[N];
    for (int i = 0; i < N; ++i) {
        this->grid_index[i] = new int[8];
        for (int j = 0; j < 8; ++j) {
            this->grid_index[i][j] = -1;
        }
    }

    this->INTERACTION_RANGE = (float) max(max(this->cohesion.sight_distance, this->alignment.sight_distance), this->separation.sight_distance);
    this->GRID_NUM = int(field_size_local_x_with_margin / (INTERACTION_RANGE*2));
    this->grid.resize(GRID_NUM*GRID_NUM*GRID_NUM);
    this->GRID_SIZE = float(field_size_local_x_with_margin / this->GRID_NUM);


    if (INTERACTION_RANGE * 2 >= this->GRID_SIZE) {
        std::cerr << "Invalid GRID_NUM!!! : " << this->GRID_NUM  << std::endl;
        exit(-1);
    }
    if (is_master) {

        cout << "# Grid List Optimization enabled" << endl
             << "# grid num :" << this->GRID_NUM << endl
             << "# grid size:" << this->GRID_SIZE << endl;
    }

    this->clear_list_grid();
}

void BoidSimulationMultiNode::clear_list_grid() {
    for (int i = 0; i < grid.size(); ++i) {
        grid[i].clear();
    }
}

void BoidSimulationMultiNode::update_list_grid() {
    this->clear_list_grid();
    for (int i = 0; i < margin_data_buffer_count/6; ++i) {
        this->update_list_grid(i);
    }
}

void BoidSimulationMultiNode::update_list_grid(unsigned int buffer_idx) {


    float local_pos[3];
    local_pos[0] = margin_data_buffer[buffer_idx*6+0]-space_x_lower+margin_width;
    local_pos[1] = margin_data_buffer[buffer_idx*6+1]-space_y_lower+margin_width;
    local_pos[2] = margin_data_buffer[buffer_idx*6+2]-space_z_lower+margin_width;
    for (int i = 0; i < 3; ++i) {
        if (local_pos[i] > field_size_local_x_with_margin) {
            local_pos[i] -= field_size;
        } else if (local_pos[i] < 0.0) {
            local_pos[i] += field_size;
        }
    }

    /*
    int grid_index_tmp[3];
    int grid_x, grid_y, grid_z, n;
    grid_index_tmp[0] = int(local_pos[0] / this->GRID_SIZE);
    grid_index_tmp[1] = int(local_pos[1] / this->GRID_SIZE);
    grid_index_tmp[2] = int(local_pos[2] / this->GRID_SIZE);
    this->grid_index[buffer_idx] = grid_index_tmp[0] + grid_index_tmp[1] * this->GRID_NUM + grid_index_tmp[2] * this->GRID_NUM * this->GRID_NUM;

    for (int i = -1; i <2; ++i) {
        grid_x = (grid_index_tmp[0] + i + GRID_NUM) % GRID_NUM;
        for (int j = -1; j <2; ++j) {
            grid_y = (grid_index_tmp[1] + j + GRID_NUM) % GRID_NUM;
            for (int k = -1; k <2; ++k) {
                grid_z = (grid_index_tmp[2] + k + GRID_NUM) % GRID_NUM;
                grid[grid_x + grid_y * GRID_NUM + grid_z * GRID_NUM * GRID_NUM].push_back(buffer_idx);
            }
        }
    }
    */

    int grid_index_tmp[3][2];
    grid_index_tmp[0][0] = int(local_pos[0] / this->GRID_SIZE);
    grid_index_tmp[1][0] = int(local_pos[1] / this->GRID_SIZE);
    grid_index_tmp[2][0] = int(local_pos[2] / this->GRID_SIZE);

    grid_index_tmp[0][0] = min(max(grid_index_tmp[0][0], 0), int(this->GRID_NUM-1));
    grid_index_tmp[1][0] = min(max(grid_index_tmp[1][0], 0), int(this->GRID_NUM-1));
    grid_index_tmp[2][0] = min(max(grid_index_tmp[2][0], 0), int(this->GRID_NUM-1));


    for (int i = 0; i < 3; ++i) {
        if (grid_index_tmp[i][0] < 0 || grid_index_tmp[i][0] >= GRID_NUM ) {


            cerr << "----debug----" << endl;
            cerr << "error2:" << buffer_idx << ":" << i << ":" << grid_index_tmp[i][0] << endl;
            cerr << "id:" << margin_data_id_buffer[buffer_idx] << endl;
            cerr << "x: " << margin_data_buffer[buffer_idx*6+0] << endl;
            cerr << "y: " << margin_data_buffer[buffer_idx*6+1] << endl;
            cerr << "z: " << margin_data_buffer[buffer_idx*6+2] << endl;
            cerr << "margin_width:" << margin_width << endl;
            cerr << "x_l:" << space_x_lower << endl;
            cerr << "y_l:" << space_y_lower << endl;
            cerr << "z_l:" << space_z_lower << endl;
            cerr << "lp_x:" << local_pos[0] << endl;
            cerr << "lp_y:" << local_pos[1] << endl;
            cerr << "lp_z:" << local_pos[2] << endl;
            cerr << "lp_z:" << local_pos[2] << endl;
            cerr << "GRID_SIZE:" << GRID_SIZE << endl;
            cerr << "rank:" << mpi_rank << endl;
            cerr << "buffer_idx:" << buffer_idx << endl;
            cerr << "field_size_local_x_with_margin:" << field_size_local_x_with_margin << endl;
            cerr << "-------------" << endl;

        }
    }


    float pos_in_grid[3];
    pos_in_grid[0] = local_pos[0] - grid_index_tmp[0][0] * this->GRID_SIZE;
    pos_in_grid[1] = local_pos[1] - grid_index_tmp[1][0] * this->GRID_SIZE;
    pos_in_grid[2] = local_pos[2] - grid_index_tmp[2][0] * this->GRID_SIZE;

    for (int i = 0; i <3; ++i) {
        if (pos_in_grid[i] <= INTERACTION_RANGE) {
            grid_index_tmp[i][1] = (grid_index_tmp[i][0] - 1 + GRID_NUM) % GRID_NUM;
        } else if (pos_in_grid[i] >= this->GRID_SIZE - INTERACTION_RANGE) {
            grid_index_tmp[i][1] = (grid_index_tmp[i][0] + 1) % GRID_NUM;
        } else {
            grid_index_tmp[i][1] = -1;
        }
    }


    // *** boid.grid_indes ***
    // first index is primally grid
    // other indices are overlapped grid and -1 mean no overlapping
    unsigned int n = 0;
    for (int i = 0; i < 2; ++i) {
        if (grid_index_tmp[0][i] == -1) continue;
        for (int j = 0; j < 2; ++j) {
            if (grid_index_tmp[1][j] == -1) continue;
            for (int k = 0; k < 2; ++k) {
                if (grid_index_tmp[2][k] == -1) continue;
                this->grid_index[buffer_idx][n] = grid_index_tmp[0][i] + grid_index_tmp[1][j] * this->GRID_NUM + grid_index_tmp[2][k] * this->GRID_NUM * this->GRID_NUM;
                n++;
            }
        }
    }


    while (n < 8) {
        this->grid_index[buffer_idx][n] = -1;
        n++;
    }

    for (int i = 0; i < 8; ++i) {
        if (this->grid_index[buffer_idx][i] < 0) {
            continue;
        }
        this->grid[this->grid_index[buffer_idx][i]].push_back(buffer_idx);
    }
}



