# MassiveSwarm

calculation of Boid model on single computer or multi node cluster (using MPI).

## Build

Mac OS X
```
$ cd ${PROJECT_ROOT}
$ mkdir build
$ cd build
$ cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ ..
```

K
```
$ cd ${PROJECT_ROOT}
$ mkdir build
$ cd build
$ cmake -D CMAKE_C_COMPILER=mpifccpx -D CMAKE_CXX_COMPILER=mpiFCCpx ..
```

## Execute

Setting file is norma ini file. Example is setting/examle.ini

```
$ boidsim -s ${setting_file} -o{output_filename} [options]
or
$ boidsim_mpi -s ${setting_file} -o{output_filename} [options]
```

### options
-h [ --help ]            help
-s [ --setting ] arg     simulation setting file with .ini format.
-o [ --output ] arg      output file.

simulation parameters (have priority over setting file.):
-N [ --population ] arg  Number of boids.
-T [ --timestep ] arg    Simulation time step.
-F [ --field-size ] arg  Size of simulation area.



## Tools

### generate job script for K

require python3

example
```
$ cd ${PROJECT_ROOT}
$ ./tools/gen_jobscript.py -s ./settings/one_body.ini ./settings/mototake.ini -n '[2**i for i in range(10,20)]' -d 16384 -t 3000 -o job_scripts/20170301
```
