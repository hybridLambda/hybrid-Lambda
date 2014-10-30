COMPILER = g++
VERSION = beta-0.3
FLAGS = -O3 -g -std=c++0x -DNDEBUG -DVERSION=\"${VERSION}\" -I../src/ 

all: clean_all hybrid-Lambda 

hybrid-Lambda: figure.o freq.o sim_gt.o
	${COMPILER} ${FLAGS} -o hybrid-Lambda debug/node_debug.cpp debug/net_debug.cpp main.cpp  mtrand.cpp  net.cpp  node.cpp hybridLambda.cpp figure.o freq.o sim_gt.o regular_math.o

sim_gt.o: regular_math.o
	${COMPILER} ${FLAGS} -c sim_gt.cpp 

regular_math.o:
	${COMPILER} ${FLAGS} -c regular_math.cpp 

freq.o:
	${COMPILER} ${FLAGS} -c freq/freq.cpp 

figure.o:
	${COMPILER} ${FLAGS} -c plot/figure.cpp

clean:
	rm -f *.o 

clean_all:
	rm -f *.o hybrid-Lambda 
