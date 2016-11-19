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
//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/ini_parser.hpp>

using std::string;
using std::cout;
using std::endl;
using std::flush;
//using boost::property_tree::ptree;
//using boost::property_tree::read_ini;

unsigned int FPS = 30;
unsigned int N;
unsigned int T;
double FIELD_SIZE;

bool is_little_endian;
BoidSimulation boid_sim;
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
         << "field size: " << boid_sim.field_size << endl
         << "#Separation" << endl
         << "force: " << boid_sim.separation.force_coefficient << endl
         << "area distance: " << boid_sim.separation.sight_distance << endl
         << "area angle: " << boid_sim.separation.sight_agnle << endl
         << "#Alignment" << endl
         << "force: " << boid_sim.alignment.force_coefficient << endl
         << "area distance: " <<  boid_sim.alignment.sight_distance << endl
         << "area angle: " << boid_sim.alignment.sight_agnle << endl
         << "#Cohesion" << endl
         << "force: " << boid_sim.cohesion.force_coefficient << endl
         << "area distance: " << boid_sim.cohesion.sight_distance << endl
         << "area angle: " << boid_sim.cohesion.sight_agnle << endl
         << "#Velocity" << endl
         << "min: " << boid_sim.velocity.min << endl
         << "max: " << boid_sim.velocity.max << endl
         << "#OpenMP: ";
    if (boid_sim.is_openmp_enabled()) {
        cout << "Enabled (max threads = " << boid_sim.get_max_threads() << ")" << endl;
    } else {
        cout << "Disabled" << endl;
    }
}

int main(int argc, char *argv[])
{
    string fname = argv[1];
    N = (unsigned int)(atoi(argv[2]));
    T = (unsigned int)(atoi(argv[3]));
    string setting_fname = argv[4];

    /*
    ptree pt;
    read_ini(setting_fname, pt);
    boid_sim.setup(N,
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
                   pt.get<double>("Velocity.MIN"));
                   */
    boid_sim.setup(N,
                   atof(argv[5]),
                   atof(argv[6]),
                   atof(argv[7]),
                   atof(argv[8]),
                   atof(argv[9]),
                   atof(argv[10]),
                   atof(argv[11]),
                   atof(argv[12]),
                   atof(argv[13]),
                   atof(argv[14]),
                   atof(argv[15]),
                   atof(argv[16]));

    print_settings();
    check_endianness();

    fout.open(fname.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    if (!fout) {
        std::cerr << "Couldn't open " << fname << endl;
        return -1;
    }

    data_file_header_v01 header;
    header.N = fix_byte_order(N);
    header.T = fix_byte_order(T);
    header.fps = fix_byte_order(FPS);
    fout.write((char*)&header,sizeof(header));

    std::cout << "simulation start." << std::endl;
    for (unsigned int t=0; t<T; t++) {
        for(int i=0; i<N; i++){
            float x = fix_byte_order((float)boid_sim.boids[i].position.x);
            float y = fix_byte_order((float)boid_sim.boids[i].position.y);
            float z = fix_byte_order((float)boid_sim.boids[i].position.z);
            fout.write((char*)&x, sizeof(x));
            fout.write((char*)&y, sizeof(y));
            fout.write((char*)&z, sizeof(z));
        }
        boid_sim.update();
        cout << t << "/" << T << "\r" << flush;
    }
    std::cout << "simulation end." << std::endl;
    fout.close();
    return 0;
}

