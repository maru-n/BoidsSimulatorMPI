# BoidsSimulatorMPI

Boids model simulator on single computer or multi node cluster (using MPI).

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
# normal version
$ boidsim -s ${setting_file} -o{output_filename} [options]

# grid optimized version
$ boidsim_grid -s ${setting_file} -o{output_filename} [options]

# MPI version
$ boidsim_mpi -s ${setting_file} -o{output_filename} [options]
```

example

- N = 16384
- time step 8000

```
./boidsim_grid -s ../settings/marun2019.ini -o output.ptcl -T 8000 -N 16384 -F 1.0
```


### options

-h [ --help ]            help

-s [ --setting ] arg     simulation setting file with .ini format.

-o [ --output ] arg      output file.

simulation parameters (have priority over setting file.):

-N [ --population ] arg  Number of boids.

-T [ --timestep ] arg    Simulation time step.

-F [ --field-size ] arg  Size of simulation area.

--vmin

--vmax

velocity limits


## Read Datafile on Python

use swarm_util.py on tools directory.

```python
import swarm_util as su
X, V = su.load_data_legacy(
    filename,
    time,
    vel_calc_step=1,
    end_time=None,
    time_interval=None,
)
```

## Tools

### generate job script for K

require python3

example
```
$ cd ${PROJECT_ROOT}
$ ./tools/gen_jobscript.py -s ./settings/one_body.ini ./settings/mototake.ini -n '[2**i for i in range(10,20)]' -d 16384 -t 3000 -o job_scripts/20170301
```
