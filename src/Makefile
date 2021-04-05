#Root libraries
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#Variables and options
OPT           = -O2
CXX           = g++
CXXFLAGS      = -Wall -fopenmp -pthread $(OPT) $(ARCH) -fPIC -I$(ROOTSYS)/include -I/home/davidg/pyrene_toymc/src
LD            = g++
LDFLAGS       = $(OPT) -fopenmp -pthread $(ARCH)
LIBS         	= $(ROOTGLIBS) -L/usr/X11R6/lib -lX11 -lXpm -lSpectrum -lm -lMathMore
DTREE_LIBS     = -L /home/davidg/pyrene_toymc/src -ltoymc

#include all header files for linking here
HDRS	= InverseCDF.h ToyMC.h

#Objects files for each class
OBJS	=	$(HDRS:.h=.o)
#Source files for each class
SRCS	=	$(HDRS:.h=.cpp)

.PHONY: all clean

all: lib Dict.C toyMC

lib: libtoymc.so

toyMC: main.o
	@echo "Creating executable.."
	$(LD) $(LDFLAGS) -o $@ $< $(LIBS) $(DTREE_LIBS)

libtoymc.so: $(OBJS) Dict.o
	$(LD) $(LDFLAGS) -shared -o libtoymc.so $^ $(LIBS)

Dict.o: Dict.C
	$(CXX) -c $(CXXFLAGS) $<

Dict.C: $(HDRS) LinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CXXFLAGS) -p $^

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $<

clean:
	rm -f *.o
	rm -f *~
	rm -f *.so
	rm -f mDict.*
	rm -f toyMC
