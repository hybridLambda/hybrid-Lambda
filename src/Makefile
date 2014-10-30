
COMPILER = g++
FLAGS = -O3 -g -std=c++0x -DNDEBUG

all: clean_all hybrid-Lambda clean

hybrid-Lambda:
	${COMPILER} ${FLAGS} -o hybrid-Lambda debug/node_debug.cpp debug/net_debug.cpp freq/freq.cpp plot/figure.cpp regular_math.cpp hybridLambda.cpp  main.cpp  mtrand.cpp  net.cpp  node.cpp  sim_gt.cpp 

clean:
	rm -f *.o 

clean_all:
	rm -f *.o hybrid-Lambda 
