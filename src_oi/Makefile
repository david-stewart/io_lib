# For now, used to compile 
# The libraries for pythia, fastjet, and root ( "LIB_TRI" )
# ccflg=-I${PYTHIA8}/include `${FASTJET3}/fastjet-config --cxxflags` `root-config --cflags`
# ldflg=-L${PYTHIA8}/lib -lpythia8 `${FASTJET3}/fastjet-config --libs` `root-config --glibs`
 
# v- ccflg= `root-config --cflags` these are the -I include paths -v
io_setter=${HOME}/root_macros/io_lib/iolib-config
ccflg=`${FASTJET3}/fastjet-config --cxxflags` `root-config --cflags` `${io_setter} -I`
#                                
# v- ccflg= `root-config --cflags` these are the -L paths
ldflg=`${FASTJET3}/fastjet-config --libs` `root-config --glibs` `${io_setter} -L`

# LIB_PYTH=-I${PYTHIA8}/include -L${PYTHIA8}/lib -lpythia8
# LIB_FASTJET=`${FASTJET3}/fastjet-config --cxxflags --libs`
# LIB_ROOT=`root-config --cflags --glibs`
# LIB_DUO= ${LIB_ROOT} ${LIB_FASTJET}
# LIB_TRI= ${LIB_ROOT} ${LIB_FASTJET} `${io_setter} -L`

# compilation option
CC=g++
CFLAGS=-std=c++11 -O3 -Wno-deprecated
# CFLAGS_CHECK=-std=c++11 -O0 -Wno-deprecated -g

lib/liboiJetMaker.so: obj/oiJetMaker.o
	${CC} -shared -o lib/liboiJetMaker.so obj/oiJetMaker.o

obj/oiJetMaker.o: src_oi/oiJetMaker.cc src_oi/oiJetMaker.h
	${CC} ${ccflg}  -fPIC -c $< -o $@

# bin/main: obj/events.o \
#           obj/main.o   \
#           obj/MemTimeProgression.o \
#           obj/test_loop.o \
#           obj/trigger_count.o \
#           obj/vz_study.o \
#           obj/tower_qa_study.o
# 	${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI} 
