//
// Created by Maruyama Norihiro on 2016/12/09.
//

#include "args.h"
#include <stdlib.h>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using std::string;
using namespace boost::property_tree;
using namespace boost::program_options;


Args::Args(int argc, const char **argv)
{
    string a;

    options_description options("");
    options.add_options()
            ("help,h", "help")
            ("setting,s", value<string>(&setting_filename), "simulation setting file with .ini format.")
            ("output,o", value<string>(&output_filename), "output file.")
            ;
    options_description sim_options("simulation parameters (have priority over setting file.)");
    sim_options.add_options()
            ("initialization,i", value<string>(&output_filename), "initialization type or file name. if file name is given, simulation will start with last conditon of the file. ")
            ("population,N", value<unsigned int>(&population), "Number of boids.")
            ("time-step,T", value<unsigned int>(&time_step), "Simulation time step.")
            ("field-size,F", value<double>(&field_size), "Size of simulation area.")
            ;
    options.add(sim_options);

    variables_map values;
    try {
        store(parse_command_line(argc, argv, options), values);
        notify(values);
    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        exit(-1);
    }

    // --help, -h
    if (values.count("help")) {
        std::cout << options << std::endl;
        exit(0);
    }

    // --output -o
    if(values.count("output")) {
    }

    // --setting, -s
    // this settings may be overwrited by other each options
    if(values.count("setting")) {
        ptree pt;
        read_ini(setting_filename, pt);
        if (boost::optional<unsigned int> o = pt.get_optional<unsigned int>("Global.TIME_STEP")) {
            time_step = o.get();
        }
        if (boost::optional<unsigned int> o = pt.get_optional<unsigned int>("Global.POPULATION")) {
            population = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Global.FIELD_SIZE")) {
            field_size = o.get();
        }
        if (boost::optional<string> o = pt.get_optional<string>("Global.INIT")) {
            initialization = o.get();
        }
        if (boost::optional<int> o = pt.get_optional<int>("Global.RANDOM_SEED")) {
            random_seed = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Separation.SIGHT_DISTANCE")) {
            separation_sight_distance = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Separation.SIGHT_ANGLE")) {
            separation_sight_angle = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Separation.FORCE_COEFFICIENT")) {
            separation_force_coefficient = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Alignment.SIGHT_DISTANCE")) {
            alignment_sight_distance = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Alignment.SIGHT_ANGLE")) {
            alignment_sight_angle = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Alignment.FORCE_COEFFICIENT")) {
            alignment_force_coefficient = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Cohesion.SIGHT_DISTANCE")) {
            cohesion_sight_distance = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Cohesion.SIGHT_ANGLE")) {
            cohesion_sight_angle = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Cohesion.FORCE_COEFFICIENT")) {
            cohesion_force_coefficient = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Velocity.MAX")) {
            velocity_max = o.get();
        }
        if (boost::optional<double> o = pt.get_optional<double>("Velocity.MIN")) {
            velocity_min = o.get();
        }
    }
    if (values.count("initialization")) {
        initialization = values["initialization"].as<string>();
    }
    if (values.count("field-size")) {
        field_size = values["field-size"].as<double>();
    }

    if (values.count("population")) {
        population = values["population"].as<unsigned int>();
    }

    if (values.count("time-step")) {
        time_step = values["time-step"].as<unsigned int>();
    }
}
