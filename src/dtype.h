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

#endif //MASSIVESWARM_DTYPE_H
