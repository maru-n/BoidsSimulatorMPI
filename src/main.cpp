/*
 *  Created by maruyama on 16/11/10.
   */

#include <stdlib.h>
#include <list>
#include <math.h>
#include <fstream>
#include <iostream>
#include "dtype.h"
#include "boid_simulation.h"
#include "boid_simulation_mpi.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using std::string;
using std::cout;
using std::endl;
using std::flush;
using boost::property_tree::ptree;
using boost::property_tree::read_ini;

unsigned int FPS = 30;
unsigned int N;
unsigned int T;
double FIELD_SIZE;

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
    if( *(char *)&x ){
        is_little_endian = true;
    }else{
        is_little_endian = false;
    }
}


void print_settings()
{
    cout << "N: " << N << endl
         << "T: " << T << endl
         << "field size: " << boid_sim->field_size_X << "," << boid_sim->field_size_Y << "," << boid_sim->field_size_Z << endl
         << "#Separation" << endl
         << "force: " << boid_sim->separation.force_coefficient << endl
         << "area distance: " << boid_sim->separation.sight_distance << endl
         << "area angle: " << boid_sim->separation.sight_agnle << endl
         << "#Alignment" << endl
         << "force: " << boid_sim->alignment.force_coefficient << endl
         << "area distance: " <<  boid_sim->alignment.sight_distance << endl
         << "area angle: " << boid_sim->alignment.sight_agnle << endl
         << "#Cohesion" << endl
         << "force: " << boid_sim->cohesion.force_coefficient << endl
         << "area distance: " << boid_sim->cohesion.sight_distance << endl
         << "area angle: " << boid_sim->cohesion.sight_agnle << endl
         << "#Velocity" << endl
         << "min: " << boid_sim->velocity.min << endl
         << "max: " << boid_sim->velocity.max << endl
         << "#OpenMP: ";
    if (boid_sim->is_openmp_enabled()) {
        cout << "Enabled (max threads = " << boid_sim->get_max_threads() << ")" << endl;
    } else {
        cout << "Disabled" << endl;
    }
}


void setup_boid_simulation(BoidSimulation* boid_sim, int argc, char *argv[])
{
    N = (unsigned int)(atoi(argv[2]));
    T = (unsigned int)(atoi(argv[3]));
    string setting_fname = argv[4];

    ptree pt;
    read_ini(setting_fname, pt);
    boid_sim->setup(N,
                    pt.get<double>("Global.FIELD_SIZE"),
                    pt.get<double>("Separation.SIGHT_DISTANCE"),
                    pt.get<double>("Separation.SIGHT_ANGLE"),
                    pt.get<double>("Separation.FORCE_COEFFICIENT"),
                    pt.get<double>("Alignment.SIGHT_DISTANCE"),
                    pt.get<double>("Alignment.SIGHT_ANGLE"),
                    pt.get<double>("Alignment.FORCE_COEFFICIENT"),
                    pt.get<double>("Cohesion.SIGHT_DISTANCE"),
                    pt.get<double>("Cohesion.SIGHT_ANGLE"),
                    pt.get<double>("Cohesion.FORCE_COEFFICIENT"),
                    pt.get<double>("Velocity.MAX"),
                    pt.get<double>("Velocity.MIN"),
                    pt.get<string>("Global.INIT"),
                    pt.get<int>("Global.RANDOM_SEED"));
}


void setup_output_data_file(string fname)
{
    check_endianness();
    fout.open(fname.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    if (!fout) {
        std::cerr << "Couldn't open output file." << endl;
        exit(-1);
    }
    data_file_header_v01 header;
    header.N = fix_byte_order(N);
    header.T = fix_byte_order(T);
    header.fps = fix_byte_order(FPS);
    fout.write((char*)&header,sizeof(header));
}


int main_singlenode(int argc, char *argv[])
{
    boid_sim = new BoidSimulation();

    setup_boid_simulation(boid_sim, argc, argv);
    print_settings();
    setup_output_data_file(argv[1]);

    // simulation
    boid_sim->init();
    std::cout << "simulation start." << std::endl;
    double x, y, z;
    float fx, fy, fz;
    for (unsigned int t=0; t<T; t++) {
        for(int i=0; i<N; i++){
            boid_sim->get(i, &x, &y, &z);
            fx = fix_byte_order(float(x));
            fy = fix_byte_order(float(y));
            fz = fix_byte_order(float(z));
            fout.write((char*)&fx, sizeof(fx));
            fout.write((char*)&fy, sizeof(fy));
            fout.write((char*)&fz, sizeof(fz));
        }
        boid_sim->update();
        cout << " " << t << "/" << T << "\r" << flush;
    }
    std::cout << "simulation end." << std::endl;
    fout.close();
    delete boid_sim;
    return 0;
}

#ifdef _MPI
int main_multinode(int argc, char *argv[])
{
    boid_sim = new BoidSimulationMultinode(argc, argv);
    setup_boid_simulation(boid_sim, argc, argv);
    if(dynamic_cast<BoidSimulationMultinode*>(boid_sim)->is_master_node()) {
        print_settings();
        setup_output_data_file(argv[1]);
    }

    // simulation
    boid_sim->init();
    if(dynamic_cast<BoidSimulationMultinode*>(boid_sim)->is_master_node()) {
        std::cout << "simulation start." << std::endl;
    }
    double x, y, z;
    float fx, fy, fz;
    for (unsigned int t=0; t<T; t++) {
        dynamic_cast<BoidSimulationMultinode*>(boid_sim)->gather_data();
        if(dynamic_cast<BoidSimulationMultinode*>(boid_sim)->is_master_node()) {
            for (int i = 0; i < N; i++) {
                boid_sim->get(i, &x, &y, &z);
                fx = fix_byte_order(float(x));
                fy = fix_byte_order(float(y));
                fz = fix_byte_order(float(z));
                fout.write((char*)&fx, sizeof(fx));
                fout.write((char*)&fy, sizeof(fy));
                fout.write((char*)&fz, sizeof(fz));
            }
            //cout << t << "/" << T << "\r" << flush;
            cout << t << "/" << T << endl;
        }
        boid_sim->update();
    }
    if(dynamic_cast<BoidSimulationMultinode*>(boid_sim)->is_master_node()) {
        std::cout << "simulation end." << std::endl;
    }
    delete boid_sim;
    return 0;
}
#endif


int main(int argc, char *argv[])
{
#ifdef _MPI
    return main_multinode(argc, argv);
#else
    return main_singlenode(argc, argv);
#endif
}
