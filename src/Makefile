COMPILER = g++
VERSION = $(shell cat ../version)
FLAGS = -O3 -g -std=c++0x -DNDEBUG -DVERSION=\"${VERSION}\" -I../src/

.PHONY: all
all: clean_all checkfile hybrid-Lambda

.PHONY: checkfile
checkfile:
	./checkfiles.sh

hybrid-Lambda: figure.o freq.o sim_gt.o
	${COMPILER} ${FLAGS} -o hybrid-Lambda debug/node_debug.cpp debug/net_debug.cpp main.cpp  fastfunc.cpp mersenne_twister.cpp random_generator.cpp net.cpp  node.cpp hybridLambda.cpp figure.o freq.o sim_gt.o regular_math.o

sim_gt.o: regular_math.o
	${COMPILER} ${FLAGS} -c sim_gt.cpp

regular_math.o:
	${COMPILER} ${FLAGS} -c regular_math.cpp

freq.o:
	${COMPILER} ${FLAGS} -c freq/freq.cpp

figure.o:
	${COMPILER} ${FLAGS} -c plot/figure.cpp

.PHONY: clean
clean:
	rm -f *.o

.PHONY: clean_all
clean_all:
	rm -f *.o hybrid-Lambda
