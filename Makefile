PROGRAM	= MassiveSwarmSimCUI
OBJS	= src/main.o src/boid.o src/boid_simulation.o src/vector3D.o

CC = fccpx
CXX = FCCpx
#CC = gcc
#CXX = g++

CFLAGS   = -V -Kfast,parallel,openmp,optmsg=2 -V -Nsrc,sta -L/volume2/home/hp160264/k03378/boost_1_62_0/stage/lib -I/volume2/home/hp160264/k03378/boost_1_62_0
CXXFLAGS = -Kfast,parallel,openmp,optmsg=2 -V -Nsrc,sta -L/volume2/home/hp160264/k03378/boost_1_62_0/stage/lib -I/volume2/home/hp160264/k03378/boost_1_62_0 -std=c++11

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(PROGRAM) $(OBJS)

