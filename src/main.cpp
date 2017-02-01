/*
 *  Created by maruyama on 16/11/10.
   */

#include <stdlib.h>
#include <list>
#include <fstream>
#include <iostream>
//#include <random>
#include "dtype.h"
#include "args.h"
#include "boid_simulation.h"
#include "boid_simulation_mpi.h"

using std::string;

const unsigned int FPS = 30;

bool is_little_endian;
BoidSimulation* boid_sim;
std::ofstream fout;


template <typename T> T fix_byte_order(T value) {
    if (is_little_endian) {
        return value;
    }
    T ret;
    int size = sizeof(T);
    for (int i = 0; i < size; ++i) {
        ((char*)&ret)[i] = ((char*)&value)[size-i-1];
    }
    return ret;
}


void check_endianness() {
    int x = 1;   // 0x00000001
    is_little_endian = *(char *)&x != 0;
}

bool is_master() {
#ifdef _MPI
    return dynamic_cast<boid_sim->ulationMultinode*>(boid_sim)->is_master_node();
#else
    return true;
#endif
}

int main(int argc, char **argv)
{
#ifdef _MPI
    boid_sim = new BoidSimulationMultinode(argc, argv);
#else
    boid_sim = new BoidSimulation();
#endif
    Args args(argc, (const char **) argv);
    boid_sim->setup(args.population, args.field_size,
                    args.separation_sight_distance, args.separation_sight_angle, args.separation_force_coefficient,
                    args.alignment_sight_distance, args.alignment_sight_angle, args.alignment_force_coefficient,
                    args.cohesion_sight_distance, args.cohesion_sight_angle, args.cohesion_force_coefficient,
                    args.velocity_max, args.velocity_min,
                    args.init_condition, args.random_seed);

    if(is_master()) {
        std::cout << args.population << std::endl;

        std::cout << "N: " << boid_sim->N << std::endl
                  << "field size: " << boid_sim->field_size_X << "," << boid_sim->field_size_Y << "," << boid_sim->field_size_Z << std::endl
                  << "#Separation" << std::endl
                  << "force: " << boid_sim->separation.force_coefficient << std::endl
                  << "area distance: " << boid_sim->separation.sight_distance << std::endl
                  << "area angle: " << boid_sim->separation.sight_agnle << std::endl
                  << "#Alignment" << std::endl
                  << "force: " << boid_sim->alignment.force_coefficient << std::endl
                  << "area distance: " <<  boid_sim->alignment.sight_distance << std::endl
                  << "area angle: " << boid_sim->alignment.sight_agnle << std::endl
                  << "#Cohesion" << std::endl
                  << "force: " << boid_sim->cohesion.force_coefficient << std::endl
                  << "area distance: " << boid_sim->cohesion.sight_distance << std::endl
                  << "area angle: " << boid_sim->cohesion.sight_agnle << std::endl
                  << "#Velocity" << std::endl
                  << "min: " << boid_sim->velocity.min << std::endl
                  << "max: " << boid_sim->velocity.max << std::endl
                  << "#OpenMP: " << (boid_sim->is_openmp_enabled() ? ("enabled (max threads:" + std::to_string(boid_sim->get_max_threads()) + ")") : "disabled")
                  << std::endl;


        check_endianness();
        fout.open(args.output_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!fout) {
            std::cerr << "Couldn't open output file." << std::endl;
            exit(-1);
        }
        data_file_header_v01 header;
        header.N = fix_byte_order(args.population);
        header.T = fix_byte_order(args.time_step);
        header.fps = fix_byte_order(FPS);
        fout.write((char *) &header, sizeof(header));
    }

    // simulation
    boid_sim->init();

    std::cout << "simulation start." << std::endl;
    double x, y, z;
    float fx, fy, fz;
    for (unsigned int t=0; t<args.time_step; t++) {
        if (is_master()) {
            for (int i = 0; i < boid_sim->N; i++) {
                boid_sim->get(i, &x, &y, &z);
                fx = fix_byte_order(float(x));
                fy = fix_byte_order(float(y));
                fz = fix_byte_order(float(z));
                fout.write((char *) &fx, sizeof(fx));
                fout.write((char *) &fy, sizeof(fy));
                fout.write((char *) &fz, sizeof(fz));
            }
            std::cout << "  " << t << "/" << args.time_step << "\r" << std::flush;
        }
        boid_sim->update();
    }
    if (is_master()) {
        std::cout << "simulation finished." << std::endl;
        fout.close();
    }
    delete boid_sim;

    return 0;
}
