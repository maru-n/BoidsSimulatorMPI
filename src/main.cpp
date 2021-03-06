/*
 *  Created by maruyama on 16/11/10.
   */

#include <stdlib.h>
#include <list>
#include <fstream>
#include <iostream>
#include "dtype.h"
#include "args.h"
#include "boid_simulation.h"
#include "boid_simulation_mpi.h"
#include <boost/format.hpp>

using std::string;

const unsigned int FPS = 30;

bool is_little_endian;
BoidSimulation* boid_sim;
std::ofstream fout;
bool is_output_v;
std::ofstream v_fout;


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
    return dynamic_cast<BoidSimulationMultiNode*>(boid_sim)->is_master_node();
#else
    return true;
#endif
}

int get_node_id(){
#ifdef _MPI
    return dynamic_cast<BoidSimulationMultiNode*>(boid_sim)->get_node_id();
#else
    return 0;
#endif
}

void write_header(std::ofstream &fout, Args& args) {
    data_file_header_v02 header;
    header.header_length = fix_byte_order(sizeof(header));
    header.N     = fix_byte_order(args.population);
    header.step  = fix_byte_order(args.time_step);
    header.t_0   = fix_byte_order(boid_sim->time_step);
    header.fps   = fix_byte_order(FPS);
    header.x_min = fix_byte_order(0.0f);
    header.x_max = fix_byte_order((float)boid_sim->field_size_X);
    header.y_min = fix_byte_order(0.0f);
    header.y_max = fix_byte_order((float)boid_sim->field_size_Y);
    header.z_min = fix_byte_order(0.0f);
    header.z_max = fix_byte_order((float) boid_sim->field_size_Z);
    fout.write((char *) &header, sizeof(header));
}


int main(int argc, char **argv)
{
    check_endianness();

#ifdef _MPI
    boid_sim = new BoidSimulationMultiNode(argc, argv);
#else
    boid_sim = new BoidSimulation();
#endif
    Args args(argc, (const char **) argv);
    boid_sim->setup(args.population, args.field_size,
                    args.separation_sight_distance, args.separation_sight_angle, args.separation_force_coefficient,
                    args.alignment_sight_distance, args.alignment_sight_angle, args.alignment_force_coefficient,
                    args.cohesion_sight_distance, args.cohesion_sight_angle, args.cohesion_force_coefficient,
                    args.velocity_max, args.velocity_min,
                    args.initialization, args.random_seed);
    boid_sim->init();
    if(is_master()) {

        std::cout << "N: " << boid_sim->N << std::endl
                  << "field size: "
                  << boid_sim->field_size_X << ","
                  << boid_sim->field_size_Y << ","
                  << boid_sim->field_size_Z << std::endl
                  << "#Separation" << std::endl
                  << "force: " << boid_sim->separation.force_coefficient << std::endl
                  << "area distance: " << boid_sim->separation.sight_distance << std::endl
                  << "area angle: " << boid_sim->separation.sight_agnle << std::endl
                  << "#Alignment" << std::endl
                  << "force: " << boid_sim->alignment.force_coefficient << std::endl
                  << "area distance: " << boid_sim->alignment.sight_distance << std::endl
                  << "area angle: " << boid_sim->alignment.sight_agnle << std::endl
                  << "#Cohesion" << std::endl
                  << "force: " << boid_sim->cohesion.force_coefficient << std::endl
                  << "area distance: " << boid_sim->cohesion.sight_distance << std::endl
                  << "area angle: " << boid_sim->cohesion.sight_agnle << std::endl
                  << "#Velocity" << std::endl
                  << "min: " << boid_sim->velocity.min << std::endl
                  << "max: " << boid_sim->velocity.max << std::endl;
        //<< "#OpenMP: " << (boid_sim->is_openmp_enabled() ? ("enabled (max threads:" + std::to_string(boid_sim->get_max_threads()) + ")") : "disabled")
        if (boid_sim->is_openmp_enabled()) {
            std::cout << "#OpenMP: enabled (max threads:" << boid_sim->get_max_threads() << ")" << std::endl;
        } else {
            std::cout << "#OpenMP: disabled" << std::endl;
        }
    }


    // setup output file
    if (args.velocity_output_filename.empty()) {
        is_output_v = false;
    } else {
        is_output_v = true;
    }

    if (args.is_parallel_output) {
        if (is_master()) {
            string fname = args.output_filename + "__header";
            std::ofstream header_file(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            if (!header_file) {
                std::cerr << "Couldn't open: " << fname << std::endl;
                exit(-1);
            }
            write_header(header_file, args);
            header_file.close();
        }
        string fname = args.output_filename + (boost::format("__c_%06d") % get_node_id()).str();
        fout.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!fout) {
            std::cerr << "Couldn't open: " << fname << std::endl;
            exit(-1);
        }
        if (is_output_v) {
            fname = args.velocity_output_filename + (boost::format("__c_%06d") % get_node_id()).str();
            v_fout.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            if (!v_fout) {
                std::cerr << "Couldn't open: " << fname << std::endl;
                exit(-1);
            }
        }
    } else {
        string fname = args.output_filename;
        fout.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!fout) {
            std::cerr << "Couldn't open: " << fname << std::endl;
            exit(-1);
        }
        write_header(fout, args);

        if (is_output_v) {
            fname = args.velocity_output_filename;
            v_fout.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            if (!v_fout) {
                std::cerr << "Couldn't open: " << fname << std::endl;
                exit(-1);
            }
            write_header(v_fout, args);
        }
    }


    // simulation
    if(is_master()) {
        std::cout << "simulation start." << std::endl;
    }
    //double x, y, z;
    //float fx, fy, fz;
    unsigned id, n, nn;
    float x, y, z, vx, vy, vz;
    float coh_x, coh_y, coh_z, sep_x, sep_y, sep_z, ali_x, ali_y, ali_z;
    for (unsigned int t=0; t<args.time_step; t++) {

        /* partial all outputs */
        unsigned output_start = 0;
        unsigned output_interval = 1;
        unsigned output_duration = 1;
        if ( (t >= output_start) && ((t - output_start) % output_interval < output_duration) ) {

            if (args.is_parallel_output) {
#ifdef _MPI
                n = dynamic_cast<BoidSimulationMultiNode*>(boid_sim)->get_data_num();
                nn = fix_byte_order(n);
                fout.write((char *) &nn, sizeof(nn));
                if (is_output_v) {
                    v_fout.write((char *) &nn, sizeof(nn));

                }
                unsigned* idb = dynamic_cast<BoidSimulationMultiNode*>(boid_sim)->get_id_buffer_ptr();
                for (int i = 0; i < n; ++i) {
                    id = idb[i];
                    id = fix_byte_order(id);
                    fout.write((char *) &id, sizeof(id));
                    if (is_output_v) {
                        v_fout.write((char *) &id, sizeof(id));
                    }
                }
                float* d = dynamic_cast<BoidSimulationMultiNode*>(boid_sim)->get_data_buffer_ptr();
                for (int i = 0; i < n; ++i) {
                    x = d[i*6+0];
                    y = d[i*6+1];
                    z = d[i*6+2];
                    x = fix_byte_order(x);
                    y = fix_byte_order(y);
                    z = fix_byte_order(z);
                    fout.write((char *) &x, sizeof(x));
                    fout.write((char *) &y, sizeof(y));
                    fout.write((char *) &z, sizeof(z));

                    if (is_output_v) {
                        x = d[i*6+3];
                        y = d[i*6+4];
                        z = d[i*6+5];
                        x = fix_byte_order(x);
                        y = fix_byte_order(y);
                        z = fix_byte_order(z);
                        v_fout.write((char *) &x, sizeof(x));
                        v_fout.write((char *) &y, sizeof(y));
                        v_fout.write((char *) &z, sizeof(z));
                    }
                }
#else
                std::cout << "Parallel output work on only MPI." << std::endl;
                return -1;
#endif
            } else {
#ifdef _MPI
                dynamic_cast<BoidSimulationMultiNode*>(boid_sim)->gather_data();
#endif
                if (is_master()) {
                    for (int i = 0; i < boid_sim->N; i++) {
                        boid_sim->get(i, &x, &y, &z);
                        x = fix_byte_order(x);
                        y = fix_byte_order(y);
                        z = fix_byte_order(z);
                        fout.write((char *) &x, sizeof(x));
                        fout.write((char *) &y, sizeof(y));
                        fout.write((char *) &z, sizeof(z));

                        if (is_output_v) {
                            boid_sim->get_v(i, &x, &y, &z);
                            x = fix_byte_order(x);
                            y = fix_byte_order(y);
                            z = fix_byte_order(z);
                            v_fout.write((char *) &x, sizeof(x));
                            v_fout.write((char *) &y, sizeof(y));
                            v_fout.write((char *) &z, sizeof(z));
                        }

                        if (args.is_force_data_output) {
                            boid_sim->get_force(i, &coh_x, &coh_y, &coh_z, &sep_x, &sep_y, &sep_z, &ali_x, &ali_y,
                                                &ali_z);
                            fout.write((char *) &coh_x, sizeof(coh_x));
                            fout.write((char *) &coh_y, sizeof(coh_y));
                            fout.write((char *) &coh_z, sizeof(coh_z));
                            fout.write((char *) &sep_x, sizeof(sep_x));
                            fout.write((char *) &sep_y, sizeof(sep_y));
                            fout.write((char *) &sep_z, sizeof(sep_z));
                            fout.write((char *) &ali_x, sizeof(ali_x));
                            fout.write((char *) &ali_y, sizeof(ali_y));
                            fout.write((char *) &ali_z, sizeof(ali_z));
                        }
                    }
                }
            }
            fout.flush();
            if (is_output_v) {
                v_fout.flush();
            }
        }

        if(is_master()) {
            //std::cout << "  " << t << "/" << args.time_step << "\r" << std::flush;
            std::cout << "  " << t << "/" << args.time_step << std::endl << std::flush;
        }
        boid_sim->update();
    }
    if (is_master()) {
        std::cout << "simulation finished." << std::endl;
        fout.close();
        if (is_output_v) {
            v_fout.close();
        }
    }
    delete boid_sim;

    return 0;
}
