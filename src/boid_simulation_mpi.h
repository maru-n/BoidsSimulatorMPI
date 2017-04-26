//
// Created by Maruyama Norihiro on 2016/11/22.
//

#ifndef MASSIVESWARM_BOID_SIMULATION_MPI_H
#define MASSIVESWARM_BOID_SIMULATION_MPI_H

#include "boid_simulation.h"

class BoidSimulationMultiNode: public BoidSimulation {
public:
    BoidSimulationMultiNode(int argc, char **argv);
    ~BoidSimulationMultiNode();
    void init();
    void update();
    //int get(unsigned int id, double* x, double* y, double* z);
    int get(unsigned int id, float* x, float* y, float* z);
    //void set_master(bool master);
    bool is_master_node(){return is_master;};

protected:
    void gather_data();

    Boid this_boid, that_boid;

    unsigned int data_buffer_count;
    //double* data_buffer;
    //double* data_buffer_swap;
    float* data_buffer;
    float* data_buffer_swap;
    unsigned* data_id_buffer;
    unsigned* data_id_buffer_swap;

    unsigned int margin_data_buffer_count;
    //double* margin_data_buffer;
    //double* margin_data_buffer_swap;
    float* margin_data_buffer;
    float* margin_data_buffer_swap;
    unsigned* margin_data_id_buffer;
    unsigned* margin_data_id_buffer_swap;

    //double *send_data_buffer;
    //double *recv_data_buffer;
    float *send_data_buffer;
    float *recv_data_buffer;
    unsigned *send_data_id_buffer;
    unsigned *recv_data_id_buffer;

    unsigned int* data_num_buffer;

    double space_x_lower, space_x_upper, space_y_lower, space_y_upper, space_z_lower, space_z_upper;
    double margin_x_lower, margin_x_upper, margin_y_lower, margin_y_upper, margin_z_lower, margin_z_upper;
    //double padding_x_lower, padding_x_upper, padding_y_lower, padding_y_upper, padding_z_lower, padding_z_upper;
    double field_size_local_x, field_size_local_y, field_size_local_z;
    double margin_width;
    bool is_master;
    int mpi_dim;
    int mpi_rank, mpi_size;
    int mpi_position_x, mpi_position_y, mpi_position_z;
    int mpi_topology_x, mpi_topology_y, mpi_topology_z;
    int neighborhood_rank[3*3*3];
    int neighborhood_num;
};

#endif //MASSIVESWARM_BOID_SIMULATION_MPI_H
