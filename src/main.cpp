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
#include "parameter.h"

using std::string;
using std::cout;
using std::endl;
using std::flush;

bool is_little_endian;
unsigned int N;
unsigned int T;
unsigned int FPS = 30;
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
    cout << "N = " << N << endl
         << "T = " << T << endl
         << "#Separation" << endl
         << "force: " << COEFF_SEPARATION << endl
         << "area distance: " << SIGHT_DISTANCE_SEPARATION << endl
         << "area angle: " << SIGHT_ANGLE_SEPARATION << endl

         << "#Alignment" << endl
         << "force: " << COEFF_ALIGNMENT << endl
         << "area distance: " << SIGHT_DISTANCE_ALIGNMENT << endl
         << "area angle: " << SIGHT_ANGLE_ALIGNMENT << endl

         << "#Cohesion" << endl
         << "force: " << COEFF_COHESION << endl
         << "area distance: " << SIGHT_DISTANCE_COHESION << endl
         << "area angle: " << SIGHT_ANGLE_COHESION << endl

         << "#Velocity" << endl
         << "min: " << MIN_VELOCITY << endl
         << "max: " << MAX_VELOCITY << endl;
    if (boid_sim.is_openmp_enabled()) {
        cout << "#OpenMP: Enabled (max threads = " << boid_sim.get_max_threads() << ")" << endl;
    } else {
        cout << "#OpenMP: Disabled" << endl;
    }
}

int main(int argc, char *argv[])
{
    string fname = argv[1];
    N = (unsigned int)(atoi(argv[2]));
    T = (unsigned int)(atoi(argv[3]));

    boid_sim.init(N);

    print_settings();
    check_endianness();

    fout.open(fname.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    if (!fout) {
        cout << "Couldn't open " << fname << endl;
        return -1;
    }

    data_file_header_v01 header;
    header.N = fix_byte_order(N);
    header.T = fix_byte_order(T);
    header.fps = fix_byte_order(FPS);
    fout.write((char*)&header,sizeof(header));

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

    fout.close();
    std::cout << "simulation end." << std::endl;

    return 0;
}
