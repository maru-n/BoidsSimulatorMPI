//
// Created by Maruyama Norihiro on 2016/12/09.
//

#ifndef MASSIVESWARM_ARGS_H
#define MASSIVESWARM_ARGS_H

#include <iostream>
#include <string>

class Args {
public:
    Args(int argc, const char **argv);
    std::string setting_filename;
    std::string output_filename;

    unsigned int time_step;

    unsigned int population;
    double field_size;
    double separation_sight_distance;
    double separation_sight_angle;
    double separation_force_coefficient;
    double alignment_sight_distance;
    double alignment_sight_angle;
    double alignment_force_coefficient;
    double cohesion_sight_distance;
    double cohesion_sight_angle;
    double cohesion_force_coefficient;
    double velocity_max;
    double velocity_min;
    std::string init_condition;
    int random_seed;
};

#endif //MASSIVESWARM_ARGS_H
