/*
 *  Created by maruyama on 16/11/10.
   */

#include <stdlib.h>
#include <list>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <fstream>
#include "boid.h"
#include <iostream>
#include "parameter.h"

using namespace std;

static const double PI = 6*asin( 0.5 );

unsigned int N;
unsigned int T;

Boid* boids;
Vector3D *dv;
unsigned long int timeStep = 0;
//const unsigned int oneStep = 200000;
unsigned long int ind = 0;

//std::ofstream ofData;
//std::ofstream ofParam;

int  dragMouseL = 0;
int  lastMouseX;
int  lastMouseY;

bool is_little_endian;
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
        std::cout << "host:little" << std::endl;
    }else{
        is_little_endian = false;
        std::cout << "host:big" << std::endl;
    }
}

void update_boids()
{
    //Vector3D dv[N];
    for(int i=0; i<N; i++){
        Vector3D dv_coh;
        Vector3D dv_sep;
        Vector3D dv_ali;
        int neivers_num_coh = 0;
        int neivers_num_sep = 0;
        int neivers_num_ali = 0;
        for(int j=0; j<N; j++){
            Vector3D boids_j_pos_tmp = boids[j].position;

            if ((boids_j_pos_tmp.x - boids[i].position.x) > FIELD_SIZE/2) {
                boids_j_pos_tmp.x -= FIELD_SIZE;
            } else if (boids[i].position.x - boids_j_pos_tmp.x > FIELD_SIZE/2) {
                boids_j_pos_tmp.x += FIELD_SIZE;
            }
            if ((boids_j_pos_tmp.y - boids[i].position.y) > FIELD_SIZE/2) {
                boids_j_pos_tmp.y -= FIELD_SIZE;
            } else if (boids[i].position.y - boids_j_pos_tmp.y > FIELD_SIZE/2) {
                boids_j_pos_tmp.y += FIELD_SIZE;
            }
            if ((boids_j_pos_tmp.z - boids[i].position.z) > FIELD_SIZE/2) {
                boids_j_pos_tmp.z -= FIELD_SIZE;
            } else if (boids[i].position.z - boids_j_pos_tmp.z > FIELD_SIZE/2) {
                boids_j_pos_tmp.z += FIELD_SIZE;
            }

            Boid target_boid(boids_j_pos_tmp);
            if( i != j ){
                // Cohesion
                if (boids[i].isInsideCohesionArea(target_boid)){
                    neivers_num_coh ++;
                    dv_coh += target_boid.position;
                }
                // Separation
                if (boids[i].isInsideSeparationArea(target_boid)) {
                    neivers_num_sep ++;
                    dv_sep += (boids[i].position - target_boid.position).getUnity();
                }
                // Alignment
                if (boids[i].isInsideAlignmentArea(target_boid)) {
                    neivers_num_ali ++;
                    dv_ali += boids[j].velocity;
                }
            }
        }
        if (neivers_num_coh != 0) {
            dv_coh = dv_coh / neivers_num_coh - boids[i].position;
        }
        if (neivers_num_sep != 0) {
            //dv_sep = dv_sep / neivers_num_sep;
        }
        if (neivers_num_ali != 0) {
            dv_ali = dv_ali / neivers_num_ali - boids[i].velocity;
        }
        dv[i] = COEFF_COHESION*dv_coh + COEFF_SEPARATION*dv_sep +COEFF_ALIGNMENT*dv_ali;
    }

    for(int i=0; i<N; i++) {
        boids[i].velocity += dv[i];

        if(boids[i].velocity.getAbs()>0. && boids[i].velocity.getAbs()>MAX_VELOCITY){
            boids[i].velocity = boids[i].velocity.getUnity() * MAX_VELOCITY;
        }else if(boids[i].velocity.getAbs()>0. && boids[i].velocity.getAbs()<MIN_VELOCITY){
            boids[i].velocity = boids[i].velocity.getUnity() * MIN_VELOCITY;
        }

        //update boid
        boids[i].position += boids[i].velocity;

        //Boundary conditon
        if(boids[i].position.x < 0.0) {
            boids[i].position.x = FIELD_SIZE + boids[i].position.x;
        }
        if(boids[i].position.y < 0.0) {
            boids[i].position.y = FIELD_SIZE + boids[i].position.y;
        }
        if(boids[i].position.z < 0.0) {
            boids[i].position.z = FIELD_SIZE + boids[i].position.z;
        }
        if(boids[i].position.x > FIELD_SIZE) {
            boids[i].position.x = boids[i].position.x - FIELD_SIZE;
        }
        if(boids[i].position.y > FIELD_SIZE) {
            boids[i].position.y = boids[i].position.y - FIELD_SIZE;
        }
        if(boids[i].position.z > FIELD_SIZE) {
            boids[i].position.z = boids[i].position.z - FIELD_SIZE;
        }

        //output data
        /*
        if(1-LOGGING && timeStep%outputStep==0){
            ofData << ind++ << "," << i << " "
                   << boids[i].position.x * 2.0 - 1.0 << " "
                   << boids[i].position.y * 2.0 - 1.0 << " "
                   << boids[i].position.z * 2.0 - 1.0 << " ";
            if(boids[i].live){
                ofData << 0;
            }else{
                ofData << 1;
            }
            ofData << ";\n";
            if(stopPeriod && ind==10000){
                exit(0);
            }
        }
         */
    }
    /*
    if(LOGGING && timeStep%outputStep==0)
        susceptability();
        */
}

void init(void)
{
    //ofData.open(dataFileName.c_str());
    //ofParam.open(paramFileName.c_str());

    std::cout
            << "N = " << N << "\n"

            << "#Separation" << "\n"
            << "force: " << COEFF_SEPARATION << "\n"
            << "area distance: " << SIGHT_DISTANCE_SEPARATION << "\n"
            << "area angle: " << SIGHT_ANGLE_SEPARATION << "\n"

            << "#Alignment" << "\n"
            << "force: " << COEFF_ALIGNMENT << "\n"
            << "area distance: " << SIGHT_DISTANCE_ALIGNMENT << "\n"
            << "area angle: " << SIGHT_ANGLE_ALIGNMENT << "\n"

            << "#Cohesion" << "\n"
            << "force: " << COEFF_COHESION << "\n"
            << "area distance: " << SIGHT_DISTANCE_COHESION << "\n"
            << "area angle: " << SIGHT_ANGLE_COHESION << "\n"

            << "#Velocity" << "\n"
            << "min: " << MIN_VELOCITY << "\n"
            << "max: " << MAX_VELOCITY << "\n";

    srand(12345);
    for(int i=0; i<N; i++){
        boids[i].position.x =FIELD_SIZE*3/8 + drand48()*FIELD_SIZE/4;
        boids[i].position.y =FIELD_SIZE*3/8 + drand48()*FIELD_SIZE/4;
        if (drand48() > 0.5) {
            boids[i].position.z = drand48()*FIELD_SIZE/4 + FIELD_SIZE*3/4;
        } else {
            boids[i].position.z = drand48()*FIELD_SIZE/4;
        }

        double v = drand48() * (MAX_VELOCITY - MIN_VELOCITY) + MIN_VELOCITY;
        double th1 = drand48() * PI;
        double th2 = drand48() * 2.0 * PI;
        boids[i].velocity.x = v * sin(th1) * cos(th2);
        boids[i].velocity.y = v * sin(th1) * sin(th2);
        boids[i].velocity.z = v * cos(th1);
    }
}

int main(int argc, char *argv[])
{
    string fname;
    fname = argv[1];
    N = (unsigned int)(atoi(argv[2]));
    T = (unsigned int)(atoi(argv[3]));

    boids = new Boid[N];
    dv = new Vector3D[N];

    check_endianness();

    ofstream fout;
    fout.open(fname.c_str(), ios::out|ios::binary|ios::trunc);
    if (!fout) {
        cout << "Couldn't open " << fname << endl;
        return -1;
    }

    std::cout << "start" << std::endl;
    init();
    unsigned char MAJOR_VERSION = 0;
    unsigned char MINOR_VERSION = 1;
    char type[] = "PTCL";
    fout.write(type, 4);

    fout.write((char*)&MAJOR_VERSION, sizeof(MAJOR_VERSION));
    fout.write((char*)&MINOR_VERSION, sizeof(MINOR_VERSION));

    unsigned int tmp;
    tmp = fix_byte_order(N);
    fout.write((char*)&tmp,sizeof(tmp));
    tmp = fix_byte_order(T);
    fout.write((char*)&tmp,sizeof(tmp));
    tmp = fix_byte_order(fps);
    fout.write((char*)&tmp,sizeof(tmp));

    for (unsigned int t=0; t<T; t++) {
        std::cout << t << std::endl;
        for(int i=0; i<N; i++){
            float x = fix_byte_order((float)boids[i].position.x);
            float y = fix_byte_order((float)boids[i].position.y);
            float z = fix_byte_order((float)boids[i].position.z);
            fout.write((char*)&x, sizeof(x));
            fout.write((char*)&y, sizeof(y));
            fout.write((char*)&z, sizeof(z));
        }
        update_boids();
    }
    //glutMainLoop();
    fout.close();
    std::cout << "finish" << std::endl;

    return 0;
}



