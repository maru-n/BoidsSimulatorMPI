//
// Created by Maruyama Norihiro on 2016/11/22.
//

#ifndef MASSIVESWARM_BOID_SIMULATION_MPI_H
#define MASSIVESWARM_BOID_SIMULATION_MPI_H

#include "boid_simulation.h"

class BoidSimulationMultinode: public BoidSimulation {
public:
    BoidSimulationMultinode(int argc, char **argv);
    ~BoidSimulationMultinode();
    void init();
    void update();
    int get(unsigned int id, double* x, double* y, double* z);
    void set_master(bool master);
    bool is_master_node(){return is_master;};

    double* data_buffer;
    unsigned* data_id_buffer;
protected:
    void gather_data();

    double* data_buffer_swap;
    unsigned* data_id_buffer_swap;
    double* margin_data_buffer;
    double* margin_data_buffer_swap;
    unsigned int margin_data_buffer_count;
    unsigned int data_buffer_count;
    unsigned int* data_num_buffer;
    double space_x_lower, space_x_upper, space_y_lower, space_y_upper, space_z_lower, space_z_upper;
    double margin_x_lower, margin_x_upper, margin_y_lower, margin_y_upper, margin_z_lower, margin_z_upper;
    double padding_x_lower, padding_x_upper, padding_y_lower, padding_y_upper, padding_z_lower, padding_z_upper;
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
