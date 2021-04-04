#Root libraries
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#Variables and options
OPT           = -O2
CXX           = g++
CXXFLAGS      = -Wall -fopenmp -pthread $(OPT) $(ARCH) -fPIC -I$(ROOTSYS)/include
LD            = g++
LDFLAGS       = $(OPT) -fopenmp -pthread $(ARCH)
LIBS         	= $(ROOTGLIBS) -L/usr/X11R6/lib -lX11 -lXpm -lSpectrum
DTREE_LIBS     = -L /home/davidg/pyrene_toymc/src -ltree

#include all header files for linking here
HDRS	= src/ToyMC.h src/InverseCDF.h

#Objects files for each class
OBJS	=	$(HDRS:.h=.o)
#Source files for each class
SRCS	=	$(HDRS:.h=.cpp)

.PHONY: all clean

all: lib Dict.C toyMC

lib: libtree.so

toyMC: $(OBJS)
	@echo "Creating executable.."
	$(LD) $(LDFLAGS) -o $@ $< $(LIBS) $(DTREE_LIBS)

libtree.so: $(OBJS) Dict.o
	$(LD) $(LDFLAGS) -shared -o libtree.so $^ $(LIBS)

Dict.o: Dict.C
	$(CXX) -c $(CXXFLAGS) $<

Dict.C: $(HDRS) LinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CXXFLAGS) -p $^

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $<

clean:
	rm -f src/*.o
	rm -f src/*~
	rm -f src/*.so
	rm -f src/Dict.*
	rm -f toyMC
