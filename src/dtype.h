//
// Created by Maruyama Norihiro on 2016/11/20.
//

#ifndef MASSIVESWARM_DTYPE_H
#define MASSIVESWARM_DTYPE_H

#pragma pack(1)
typedef struct {
    char format_name[4] = {'P', 'T', 'C', 'L'};
    unsigned char major_version = 0;
    unsigned char minor_version = 1;
    unsigned int N;
    unsigned int T;
    unsigned int fps;
} data_file_header_v01;
#pragma pack()

#pragma pack(1)
typedef struct {
    char format_name[4] = {'P', 'T', 'C', 'L'};
    unsigned char major_version = 0;
    unsigned char minor_version = 2;
    unsigned int header_length;
    unsigned int N;
    unsigned int step;
    unsigned int t_0;
    unsigned int fps;
    float x_min;
    float x_max;
    float y_min;
    float y_max;
    float z_min;
    float z_max;
} data_file_header_v02;
#pragma pack()


typedef struct {
    double sight_distance;
    double sight_agnle;
    double force_coefficient;
} interaction_parameters;


typedef struct {
    double max;
    double min;
} velocity_parameters;





#endif //MASSIVESWARM_DTYPE_H
