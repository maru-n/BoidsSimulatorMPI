//
// Created by Maruyama Norihiro on 2016/11/22.
//

#include "boid.h"
#include "boid_simulation_mpi.h"
#include "mpi.h"
#include <mpi-ext.h>

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
    /*
    std::cout << "rank:" << mpi_rank << "(" << mpi_position_x << "," << mpi_position_y << "," << mpi_position_z << ")";
    for (int i = 0; i < neighborhood_num; ++i) {
        std::cout << neighborhood_rank[i] << ",";
    }
    std::cout << std::endl;
     */
}

BoidSimulationMultiNode::~BoidSimulationMultiNode()
{
    MPI_Finalize();
    delete[] data_num_buffer;
    delete[] data_buffer;
    delete[] data_id_buffer;
    delete[] data_id_buffer_swap;
    delete[] data_buffer_swap;
    delete[] margin_data_buffer;
    delete[] margin_data_buffer_swap;
}

void BoidSimulationMultiNode::init()
{
    margin_data_buffer = new double[N*6];
    margin_data_buffer_swap = new double[N*6];
    data_buffer = new double[N*6];
    data_id_buffer = new unsigned[N];
    data_buffer_swap = new double[N*6];
    data_id_buffer_swap = new unsigned[N];
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

    unsigned id;
    double x, y, z, vx, vy, vz;
    data_buffer_count = 0;
    margin_data_buffer_count = 0;
    for (int i = 0; i < N; ++i) {
        id = i;
        x = data_buffer[i*6+0];
        y = data_buffer[i*6+1];
        z = data_buffer[i*6+2];
        vx = data_buffer[i*6+3];
        vy = data_buffer[i*6+4];
        vz = data_buffer[i*6+5];
        boids[i].set(x, y, z, vx, vy, vz);

        if ( (margin_x_lower < margin_x_upper ? (x>=margin_x_lower&&x<=margin_x_upper) : (x>=margin_x_lower||x<=margin_x_upper)) &&
             (margin_y_lower < margin_y_upper ? (y>=margin_y_lower&&y<=margin_y_upper) : (y>=margin_y_lower||y<=margin_y_upper)) &&
             (margin_z_lower < margin_z_upper ? (z>=margin_z_lower&&z<=margin_z_upper) : (z>=margin_z_lower||z<=margin_z_upper))) {
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
            data_id_buffer[data_buffer_count/6] = id;
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

    this->gather_data();
}


void BoidSimulationMultiNode::update()
{
    /*
     * Update
     */

    double x, y, z, vx, vy, vz;


#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp parallel
#endif

#if defined(_OPENMP) && defined(ENABLE_OPENMP)
#pragma omp for
//#pragma omp for schedule(guided)
#endif
    for(int i=0; i<data_buffer_count/6; i++) {
        dv_coh[i].x = dv_coh[i].y = dv_coh[i].z =
        dv_sep[i].x = dv_sep[i].y = dv_sep[i].z =
        dv_ali[i].x = dv_ali[i].y = dv_ali[i].z = 0.0;
        int neivers_num_coh = 0;
        int neivers_num_sep = 0;
        int neivers_num_ali = 0;

        boids[i].set_serialized_data(&data_buffer[i*6]);
        for(int j=0; j<margin_data_buffer_count/6; j++) {
            boids[j].set_serialized_data(&margin_data_buffer[j*6]);


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
            boids[i].position.x = field_size + boids[i].position.x;
        } else if (boids[i].position.x >= field_size) {
            boids[i].position.x = boids[i].position.x - field_size;
        }
        if (boids[i].position.y < 0.0) {
            boids[i].position.y = field_size + boids[i].position.y;
        } else if (boids[i].position.y >= field_size) {
            boids[i].position.y = boids[i].position.y - field_size;
        }
        if (boids[i].position.z < 0.0) {
            boids[i].position.z = field_size + boids[i].position.z;
        } else if (boids[i].position.z >= field_size) {
            boids[i].position.z = boids[i].position.z - field_size;
        }

        boids[i].get_serialized_data(&data_buffer[i * 6]);
    }


    //Debug
    /*
    unsigned BUG_ID = 20126;
    for (int i = 0; i < data_buffer_count/6; ++i) {
        if (data_id_buffer[i] == BUG_ID) {
            std::cout << "rank:" << mpi_rank << std::endl;
            std::cout << data_buffer[i*6+0] << ","
                      << data_buffer[i*6+1] << ","
                      << data_buffer[i*6+2] << std::endl;
        }
    }
     */


    /*
     * Comunicate with other nodes
     */

    unsigned int send_data_buffer_count = 0;
    unsigned int recv_data_buffer_count = 0;
    double *send_data_buffer = new double[N*6];
    double *recv_data_buffer = new double[N*6];
    unsigned *send_data_id_buffer = new unsigned[N];
    unsigned *recv_data_id_buffer = new unsigned[N];
    //data_buffer_count = 0;

    unsigned id;
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
            data_buffer_swap[data_buffer_count_new+0] = x;
            data_buffer_swap[data_buffer_count_new+1] = y;
            data_buffer_swap[data_buffer_count_new+2] = z;
            data_buffer_swap[data_buffer_count_new+3] = vx;
            data_buffer_swap[data_buffer_count_new+4] = vy;
            data_buffer_swap[data_buffer_count_new+5] = vz;
            data_id_buffer_swap[data_buffer_count_new/6] = id;
            data_buffer_count_new += 6;
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

    unsigned* tmp2;
    tmp2 = data_id_buffer;
    data_id_buffer = data_id_buffer_swap;
    data_id_buffer_swap = tmp2;

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
        //std::cout << mpi_rank << "->" << neighborhood_rank[i] << ":" << send_data_buffer[0] << std::endl;
        MPI_Isend(send_data_buffer, send_data_buffer_count, MPI_DOUBLE, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Irecv(&recv_data_buffer[recv_data_buffer_count], data_num[i], MPI_DOUBLE, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Isend(send_data_id_buffer, send_data_buffer_count/6, MPI_UNSIGNED, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        MPI_Irecv(&recv_data_id_buffer[recv_data_buffer_count/6], data_num[i]/6, MPI_UNSIGNED, neighborhood_rank[i], 0, MPI_COMM_WORLD, &request2[request_n++]);
        recv_data_buffer_count += data_num[i];
    }

    MPI_Waitall(request_n, request2, status2);

    /*
    unsigned int a = 0;
    for (int i = 0; i < neighborhood_num; ++i) {
        //std::cout << mpi_rank << "<-" << neighborhood_rank[i] << ":" << recv_data_buffer[a] << std::endl;
        a += data_num[i];
    }
     */


    for (int i = 0; i < recv_data_buffer_count/6; i++) {
        id = recv_data_id_buffer[i];
        x = recv_data_buffer[i*6+0];
        y = recv_data_buffer[i*6+1];
        z = recv_data_buffer[i*6+2];
        vx = recv_data_buffer[i*6+3];
        vy = recv_data_buffer[i*6+4];
        vz = recv_data_buffer[i*6+5];

        if (x >= space_x_lower && x < space_x_upper &&
            y >= space_y_lower && y < space_y_upper &&
            z >= space_z_lower && z < space_z_upper) {
            data_id_buffer[data_buffer_count/6] = id;
            data_buffer[data_buffer_count+0] = x;
            data_buffer[data_buffer_count+1] = y;
            data_buffer[data_buffer_count+2] = z;
            data_buffer[data_buffer_count+3] = vx;
            data_buffer[data_buffer_count+4] = vy;
            data_buffer[data_buffer_count+5] = vz;
            data_buffer_count += 6;
            margin_data_buffer[margin_data_buffer_count+0] = x;
            margin_data_buffer[margin_data_buffer_count+1] = y;
            margin_data_buffer[margin_data_buffer_count+2] = z;
            margin_data_buffer[margin_data_buffer_count+3] = vx;
            margin_data_buffer[margin_data_buffer_count+4] = vy;
            margin_data_buffer[margin_data_buffer_count+5] = vz;
            margin_data_buffer_count +=6;
        }

    }

    for (int i = 0; i < recv_data_buffer_count/6; i++) {
        x = recv_data_buffer[i*6+0];
        y = recv_data_buffer[i*6+1];
        z = recv_data_buffer[i*6+2];
        vx = recv_data_buffer[i*6+3];
        vy = recv_data_buffer[i*6+4];
        vz = recv_data_buffer[i*6+5];

        if (!(x >= space_x_lower && x < space_x_upper &&
            y >= space_y_lower && y < space_y_upper &&
            z >= space_z_lower && z < space_z_upper)
                &&
            ( (margin_x_lower < margin_x_upper ? (x>=margin_x_lower&&x<margin_x_upper) : (x>=margin_x_lower||x<margin_x_upper)) &&
              (margin_y_lower < margin_y_upper ? (y>=margin_y_lower&&y<margin_y_upper) : (y>=margin_y_lower||y<margin_y_upper)) &&
              (margin_z_lower < margin_z_upper ? (z>=margin_z_lower&&z<margin_z_upper) : (z>=margin_z_lower||z<margin_z_upper)))){
            margin_data_buffer[margin_data_buffer_count+0] = x;
            margin_data_buffer[margin_data_buffer_count+1] = y;
            margin_data_buffer[margin_data_buffer_count+2] = z;
            margin_data_buffer[margin_data_buffer_count+3] = vx;
            margin_data_buffer[margin_data_buffer_count+4] = vy;
            margin_data_buffer[margin_data_buffer_count+5] = vz;
            margin_data_buffer_count +=6;
        }
    }



    //std::cout << mpi_rank << ": n:" << data_buffer_count << ": mn" << margin_data_buffer_count << std::endl;
    this->gather_data();

    delete[] send_data_buffer;
    delete[] recv_data_buffer;
    delete[] send_data_id_buffer;
    delete[] recv_data_id_buffer;
}

//int BoidSimulationMultiNode::get(unsigned int id, double* x, double* y, double* z)
int BoidSimulationMultiNode::get(unsigned int id, float* x, float* y, float* z)
{
    //std::cout << id << std::endl;
    for (int i = 0; i < N; ++i) {
        if(data_id_buffer[i] == id) {
            *x = data_buffer[6*i+0];
            *y = data_buffer[6*i+1];
            *z = data_buffer[6*i+2];
            return 0;
        }
    }
    std::cerr << "exit(-1) called with id:" << id << std::endl;
    exit(-1);
    /*
    *x = data_buffer[6*id+0];
    *y = data_buffer[6*id+1];
    *z = data_buffer[6*id+2];
    return 0;
     */
}

void BoidSimulationMultiNode::gather_data()
{
    MPI_Gather(&data_buffer_count, 1, MPI_INT, data_num_buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (is_master) {
        //MPI_Request request[(mpi_size-1)*2];
        //MPI_Status status[(mpi_size-1)*2];
        MPI_Request *request = new MPI_Request[(mpi_size-1)*2];
        MPI_Status *status = new MPI_Status[(mpi_size-1)*2];
        int n = data_num_buffer[0];
        //std::cerr << "data_num" << std::endl;
        //std::cerr << data_num_buffer[0]/6 << ",";
        for (int i = 1; i < mpi_size; ++i) {
            //MPI_Status status;
            MPI_Irecv(&data_buffer[n], data_num_buffer[i], MPI_DOUBLE, i, (i*2), MPI_COMM_WORLD, &request[(i-1)*2]);
            MPI_Irecv(&data_id_buffer[n/6], data_num_buffer[i]/6, MPI_UNSIGNED, i, (i*2)+1, MPI_COMM_WORLD, &request[(i-1)*2+1]);
            //std::cerr << data_num_buffer[i]/6 << ",";
            n += data_num_buffer[i];
        }
        //std::cerr << std::endl;

        MPI_Waitall((mpi_size-1)*2, request, status);
        delete[] request;
        delete[] status;

        // Debug
        //std::cerr << "data_id" << std::endl;
        /*
        for (int j = 0; j < N; ++j) {
            std::cerr << data_id_buffer[j] << ",";
        }
        std::cerr << std::endl;
         */

    } else {
        MPI_Request request[2];
        MPI_Status status[2];
        MPI_Isend(data_buffer, data_buffer_count, MPI_DOUBLE, 0, mpi_rank*2, MPI_COMM_WORLD, &request[0]);
        MPI_Isend(data_id_buffer, data_buffer_count/6, MPI_UNSIGNED, 0, (mpi_rank*2)+1, MPI_COMM_WORLD, &request[1]);
        //MPI_Send(data_buffer, data_buffer_count, MPI_DOUBLE, 0, mpi_rank*2, MPI_COMM_WORLD);
        //MPI_Send(data_id_buffer, data_buffer_count/6, MPI_UNSIGNED, 0, (mpi_rank*2)+1, MPI_COMM_WORLD);
        MPI_Waitall(2, request, status);
    }
}

