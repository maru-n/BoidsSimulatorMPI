//
// Created by Maruyama Norihiro on 2016/11/22.
//

#include "boid.h"
#include "boid_simulation_multinode.h"
#include "mpi.h"
#include <mpi-ext.h>

using namespace std;

BoidSimulationMultinode::BoidSimulationMultinode(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    FJMPI_Topology_get_dimension(&mpi_dim);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    FJMPI_Topology_get_shape(&mpi_topology_x, &mpi_topology_y, &mpi_topology_z);
    FJMPI_Topology_rank2xyz(mpi_rank, &mpi_position_x, &mpi_position_y, &mpi_position_z);
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
                //FJMPI_Topology_xyz2rank(n_x, n_y, n_z, &neighborhood_rank[i][j][k]);
                FJMPI_Topology_xyz2rank(n_x, n_y, n_z, &r);
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
                //if (is_master) std::cout << "master:n_rank:" << i << "," << j <<  "," << k << "(" << neighborhood_rank[i][j][k] << ")" << std::endl;
            }
        }
    }
}

BoidSimulationMultinode::~BoidSimulationMultinode()
{
    MPI_Finalize();
    delete[] data_num_buffer;
    delete[] data_buffer;
    delete[] data_buffer_swap;
    delete[] margin_data_buffer;
    delete[] margin_data_buffer_swap;
}

void BoidSimulationMultinode::init()
{
    margin_data_buffer = new double[N*6];
    margin_data_buffer_swap = new double[N*6];
    data_buffer = new double[N*6];
    data_buffer_swap = new double[N*6];

    data_num_buffer = new unsigned int[mpi_size];
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
            boids[i].get_serialized_data(&data_buffer[i*6]);
        }
    }
    MPI_Bcast(data_buffer, N*6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double x, y, z, vx, vy, vz;
    data_buffer_count = 0;
    margin_data_buffer_count = 0;
    for (int i = 0; i < N; ++i) {
        x = data_buffer[i*6+0];
        y = data_buffer[i*6+1];
        z = data_buffer[i*6+2];
        vx = data_buffer[i*6+3];
        vy = data_buffer[i*6+4];
        vz = data_buffer[i*6+5];
        boids[i].set(x, y, z, vx, vy, vz);

        if ( (margin_x_lower < margin_x_upper ? (x>=margin_x_lower&&x<margin_x_upper) : (x>=margin_x_lower||x<margin_x_upper)) &&
             (margin_y_lower < margin_y_upper ? (y>=margin_y_lower&&x<margin_y_upper) : (y>=margin_y_lower||y<margin_y_upper)) &&
             (margin_z_lower < margin_z_upper ? (z>=margin_z_lower&&x<margin_z_upper) : (z>=margin_z_lower||z<margin_z_upper))) {
            margin_data_buffer[margin_data_buffer_count+0] = x;
            margin_data_buffer[margin_data_buffer_count+1] = y;
            margin_data_buffer[margin_data_buffer_count+2] = z;
            margin_data_buffer[margin_data_buffer_count+3] = vx;
            margin_data_buffer[margin_data_buffer_count+4] = vy;
            margin_data_buffer[margin_data_buffer_count+5] = vz;
            margin_data_buffer_count +=6;
        }
        if (x >= space_x_lower && x < space_x_upper &&
            y >= space_y_lower && y < space_y_upper &&
            z >= space_z_lower && z < space_z_upper) {
            data_buffer_swap[data_buffer_count+0] = x;
            data_buffer_swap[data_buffer_count+1] = y;
            data_buffer_swap[data_buffer_count+2] = z;
            data_buffer_swap[data_buffer_count+3] = vx;
            data_buffer_swap[data_buffer_count+4] = vy;
            data_buffer_swap[data_buffer_count+5] = vz;
            data_buffer_count += 6;
        }
    }
    double* tmp = data_buffer;
    data_buffer = data_buffer_swap;
    data_buffer_swap = tmp;
}


void BoidSimulationMultinode::update()
{
    double x, y, z, vx, vy, vz;
#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp parallel
#endif

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
//#pragma omp for
#pragma omp for schedule(guided)
#endif
    //for(int i=0; i<local_n; i++){
    //for(int i=0; i<N; i++){
    for(int i=0; i<data_buffer_count/6; i++) {
        //if (!boid_in_local_area_flag[i]) continue;
        dv_coh[i].x = dv_coh[i].y = dv_coh[i].z =
        dv_sep[i].x = dv_sep[i].y = dv_sep[i].z =
        dv_ali[i].x = dv_ali[i].y = dv_ali[i].z = 0.0;
        int neivers_num_coh = 0;
        int neivers_num_sep = 0;
        int neivers_num_ali = 0;
        //for(int j=0; j<local_n; j++){
        //for(int j=0; j<N; j++){
        boids[i].set_serialized_data(&data_buffer[i*6]);
        for(int j=0; j<margin_data_buffer_count/6; j++) {
            boids[j].set_serialized_data(&margin_data_buffer[j*6]);
            //if (!boid_in_margin_area_flag[j]) continue;

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
//#pragma omp for
#pragma omp for schedule(guided)
#endif
    //for(int i=0; i<local_n; i++) {
    for(int i=0; i<data_buffer_count/6; i++) {
        boids[i].set_serialized_data(&data_buffer[i * 6]);
        //for(int i=0; i<N; i++){
        //if (!boid_in_local_area_flag[i]) continue;

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
            boids[i].position.x = field_size + boids[i].position.x;
        } else if (boids[i].position.x > field_size) {
            boids[i].position.x = boids[i].position.x - field_size;
        }
        if (boids[i].position.y < 0.0) {
            boids[i].position.y = field_size + boids[i].position.y;
        } else if (boids[i].position.y > field_size) {
            boids[i].position.y = boids[i].position.y - field_size;
        }
        if (boids[i].position.z < 0.0) {
            boids[i].position.z = field_size + boids[i].position.z;
        } else if (boids[i].position.z > field_size) {
            boids[i].position.z = boids[i].position.z - field_size;
        }
        boids[i].get_serialized_data(&data_buffer[i * 6]);
    }

    unsigned int send_data_buffer_count = 0;
    unsigned int recv_data_buffer_count = 0;
    double *send_data_buffer = new double[N*6];
    double *recv_data_buffer = new double[N*6];
    //data_buffer_count = 0;

    unsigned int data_buffer_count_new = 0;
    unsigned int margin_data_buffer_count_new = 0;
    for(int i=0; i<data_buffer_count/6; i++) {

        x = data_buffer[i*6+0];
        y = data_buffer[i*6+1];
        z = data_buffer[i*6+2];
        vx = data_buffer[i*6+3];
        vy = data_buffer[i*6+4];
        vz = data_buffer[i*6+5];

        if (x < (space_x_lower+margin_width) || x >= (space_x_upper-margin_width) ||
            y < (space_y_lower+margin_width) || y >= (space_y_upper-margin_width) ||
            z < (space_z_lower+margin_width) || z >= (space_z_upper-margin_width)) {
            send_data_buffer[send_data_buffer_count+0] = x;
            send_data_buffer[send_data_buffer_count+1] = y;
            send_data_buffer[send_data_buffer_count+2] = z;
            send_data_buffer[send_data_buffer_count+3] = vx;
            send_data_buffer[send_data_buffer_count+4] = vy;
            send_data_buffer[send_data_buffer_count+5] = vz;
            send_data_buffer_count+=6;

            //boid_in_local_area_flag[i] = false;
            //boid_in_margin_area_flag[i] = false;
        } else {
            //boids[i].get_serialized_data(&data_buffer[data_buffer_count]);
            data_buffer_swap[data_buffer_count_new+0] = x;
            data_buffer_swap[data_buffer_count_new+1] = y;
            data_buffer_swap[data_buffer_count_new+2] = z;
            data_buffer_swap[data_buffer_count_new+3] = vx;
            data_buffer_swap[data_buffer_count_new+4] = vy;
            data_buffer_swap[data_buffer_count_new+5] = vz;
            data_buffer_count_new += 6;
            //boid_in_local_area_flag[i] = true;
            //boid_in_margin_area_flag[i] = true;
            margin_data_buffer_swap[margin_data_buffer_count_new+0] = x;
            margin_data_buffer_swap[margin_data_buffer_count_new+1] = y;
            margin_data_buffer_swap[margin_data_buffer_count_new+2] = z;
            margin_data_buffer_swap[margin_data_buffer_count_new+3] = vx;
            margin_data_buffer_swap[margin_data_buffer_count_new+4] = vy;
            margin_data_buffer_swap[margin_data_buffer_count_new+5] = vz;
            margin_data_buffer_count_new +=6;

        }
    }
    double* tmp;
    tmp = data_buffer;
    data_buffer = data_buffer_swap;
    data_buffer_swap = tmp;
    data_buffer_count = data_buffer_count_new;

    tmp = margin_data_buffer;
    margin_data_buffer = margin_data_buffer_swap;
    margin_data_buffer_swap = tmp;
    margin_data_buffer_count = margin_data_buffer_count_new;


    MPI_Request request1[27*2];
    MPI_Request request2[27*2];
    MPI_Status status[27*2];
    int request_n=0;
    /*
    int data_num[3][3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                //std::cout << mpi_position_x << "," << mpi_position_y << "," << mpi_position_z << ":send_num:" << send_data_buffer_count << std::endl;
                MPI_Isend(&send_data_buffer_count, 1, MPI_INT, neighborhood_rank[i][j][k], 0, MPI_COMM_WORLD, &request1[request_n++]);
                MPI_Irecv(&data_num[i][j][k], 1, MPI_INT, neighborhood_rank[i][j][k], 0, MPI_COMM_WORLD, &request1[request_n++]);
            }
        }
    }

    MPI_Waitall(27*2, request1, status);
     */
    int data_num[27];
    for (int i = 0; i < neighborhood_num; ++i) {
        //if (neighborhood_rank[i]==0) std::cout << mpi_rank << ":send_to_master:"  << send_data_buffer_count << std::endl;
        MPI_Isend(&send_data_buffer_count, 1, MPI_INT, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request1[request_n++]);
        MPI_Irecv(&data_num[i], 1, MPI_INT, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request1[request_n++]);
    }
    MPI_Waitall(request_n, request1, status);

    /*
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                //if (is_master) std::cout << "master_receive:" << i << "," << j <<  "," << k << "(" << neighborhood_rank[i][j][k] << "):" << data_num[i][j][k] << std::endl;
                MPI_Isend(send_data_buffer, send_data_buffer_count, MPI_DOUBLE, neighborhood_rank[i][j][k], 0, MPI_COMM_WORLD, &request2[request_n++]);
                MPI_Irecv(&recv_data_buffer[recv_data_buffer_count], data_num[i][j][k], MPI_DOUBLE, neighborhood_rank[i][j][k], 0, MPI_COMM_WORLD, &request2[request_n++]);
                recv_data_buffer_count += data_num[i][j][k];
            }
        }
    }
    MPI_Waitall(27*2, request2, status);
     */
    request_n = 0;
    for (int i = 0; i < neighborhood_num; ++i) {
        //if (is_master) std::cout << "master_receive:" << neighborhood_rank[i] << ":" << data_num[i] << std::endl;
        MPI_Isend(send_data_buffer, send_data_buffer_count, MPI_DOUBLE, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Irecv(&recv_data_buffer[recv_data_buffer_count], data_num[i], MPI_DOUBLE, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        recv_data_buffer_count += data_num[i];
    }
    MPI_Waitall(request_n, request2, status);


    for (int i = 0; i < recv_data_buffer_count/6; i++) {
        x = recv_data_buffer[i*6+0];
        y = recv_data_buffer[i*6+1];
        z = recv_data_buffer[i*6+2];
        vx = recv_data_buffer[i*6+3];
        vy = recv_data_buffer[i*6+4];
        vz = recv_data_buffer[i*6+5];

        if ( (margin_x_lower < margin_x_upper ? (x>=margin_x_lower&&x<margin_x_upper) : (x>=margin_x_lower||x<margin_x_upper)) &&
             (margin_y_lower < margin_y_upper ? (y>=margin_y_lower&&x<margin_y_upper) : (y>=margin_y_lower||y<margin_y_upper)) &&
             (margin_z_lower < margin_z_upper ? (z>=margin_z_lower&&x<margin_z_upper) : (z>=margin_z_lower||z<margin_z_upper))) {
            margin_data_buffer[margin_data_buffer_count+0] = x;
            margin_data_buffer[margin_data_buffer_count+1] = y;
            margin_data_buffer[margin_data_buffer_count+2] = z;
            margin_data_buffer[margin_data_buffer_count+3] = vx;
            margin_data_buffer[margin_data_buffer_count+4] = vy;
            margin_data_buffer[margin_data_buffer_count+5] = vz;
            margin_data_buffer_count +=6;
        }
        if (x >= space_x_lower && x < space_x_upper &&
            y >= space_y_lower && y < space_y_upper &&
            z >= space_z_lower && z < space_z_upper) {
            data_buffer[data_buffer_count+0] = x;
            data_buffer[data_buffer_count+1] = y;
            data_buffer[data_buffer_count+2] = z;
            data_buffer[data_buffer_count+3] = vx;
            data_buffer[data_buffer_count+4] = vy;
            data_buffer[data_buffer_count+5] = vz;
            data_buffer_count += 6;
        }

    }

    delete[] send_data_buffer;
    delete[] recv_data_buffer;
}

void BoidSimulationMultinode::gather_data()
{
    MPI_Gather(&data_buffer_count, 1, MPI_INT, data_num_buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (is_master) {
        int n = data_num_buffer[0];
        for (int i = 1; i < mpi_size; ++i) {
            MPI_Status status;
            MPI_Recv(&data_buffer[n], data_num_buffer[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            n += data_num_buffer[i];
        }
        /*
        for (int i=0; i<N; i++) {
            boids[i].set_serialized_data(&data_buffer[i*6]);
        }
         */
    } else {
        MPI_Send(data_buffer, data_buffer_count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}


