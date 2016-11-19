PROGRAM	= MassiveSwarmSimCUI
OBJS	= src/main.o src/boid.o src/boid_simulation.o src/vector3D.o

CC	= fccpx
CXX = FCCpx

CFLAGS   = -Kfast,parallel,openmp,optmsg=2 -V -Nsrc,sta
CXXFLAGS = -Kfast,parallel,openmp,optmsg=2 -V -Nsrc,sta -std=c++11

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(PROGRAM) $(OBJS)

