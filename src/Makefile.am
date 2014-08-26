bin_PROGRAMS = hybrid-Lambda
PROG = HYBRID-LAMBDA
check_PROGRAMS = hybrid-Lambda_dbg hybrid-Lambda_prof unit_tests
TESTS = unit_tests 

hybrid_Lambda_CXXFLAGS = -O3 -g -DNDEBUG #options for not showing debug messages.
hybrid_Lambda_dbg_CXXFLAGS =  -O3 -g 
hybrid_Lambda_prof_CXXFLAGS = -O3 -g -DNDEBUG #options for not showing debug messages.

shared_SOURCES = hybridLambda.cpp \
				 node.cpp \
				 net.cpp \
				 sim_gt.cpp \
				 freq.cpp \
				 utility.cpp  mtrand.cpp \
				 figure.cpp  

hybrid_Lambda_SOURCES = $(shared_SOURCES) main.cpp  
hybrid_Lambda_dbg_SOURCES = $(shared_SOURCES) main.cpp  
hybrid_Lambda_prof_SOURCES = $(shared_SOURCES) main.cpp  

unit_tests_SOURCES = $(shared_SOURCES) tests/test_hybridLambda.cpp tests/test_figure.cpp tests/test_runner.cpp 

unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 
unit_tests_LDADD= -lcppunit -ldl

include_HEADERS = sim_gt.hpp figure.hpp freq.hpp utility.hpp mtrand.h hybridLambda.hpp node.hpp net.hpp 


