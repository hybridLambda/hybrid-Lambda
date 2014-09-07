EXTRA_DIST = trees/*tax* old_Makefile bootstrap example/eg* 
#DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
distdir = $(PACKAGE)-$(VERSION)


bin_PROGRAMS = hybrid-Lambda
PROG = HYBRID-LAMBDA
check_PROGRAMS = hybrid-Lambda_dbg hybrid-Lambda_prof unit_tests
TESTS = unit_tests 

hybrid_Lambda_CXXFLAGS = -O3 -g -DNDEBUG #options for not showing debug messages.
hybrid_Lambda_dbg_CXXFLAGS =  -O3 -g 
hybrid_Lambda_prof_CXXFLAGS = -O3 -g -DNDEBUG #options for not showing debug messages.

shared_SOURCES = src/hybridLambda.cpp \
				 src/node.cpp \
				 src/net.cpp \
				 src/sim_gt.cpp \
				 src/freq.cpp \
				 src/utility.cpp \
				 src/mtrand.cpp \
				 src/figure.cpp  
				 
debug_SOURCES = src/debug/net_debug.cpp \
			    src/debug/node_debug.cpp

hybrid_Lambda_SOURCES = $(shared_SOURCES) src/main.cpp  
hybrid_Lambda_dbg_SOURCES = $(shared_SOURCES) $(debug_SOURCES) src/main.cpp  
hybrid_Lambda_prof_SOURCES = $(shared_SOURCES) src/main.cpp  

unit_tests_SOURCES = $(shared_SOURCES) src/tests/test_hybridLambda.cpp src/tests/test_figure.cpp src/tests/test_runner.cpp 

unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 
unit_tests_LDADD= -lcppunit -ldl

include_HEADERS = src/sim_gt.hpp src/figure.hpp src/freq.hpp src/utility.hpp src/mtrand.h src/hybridLambda.hpp src/node.hpp src/net.hpp 

