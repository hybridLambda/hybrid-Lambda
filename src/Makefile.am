bin_PROGRAMS = hybrid-Lambda
PROG = HYBRID-LAMBDA
check_PROGRAMS = hybrid-Lambda_dbg hybrid-Lambda_prof unit_tests
TESTS = unit_tests 

hybrid_Lambda_CXXFLAGS = -O3 -g -DNDEBUG #options for not showing debug messages.
hybrid_Lambda_dbg_CXXFLAGS =  -O3 -g 
hybrid_Lambda_prof_CXXFLAGS = -O3 -g -DNDEBUG #options for not showing debug messages.


shared_SOURCES =  param.hpp param.cpp \
						node.hpp node.cpp \
						net.hpp net.cpp \
						all_gene_topo.cpp  all_gene_topo.hpp \
						sim_gt.cpp  sim_gt.hpp  \
						freq.cpp  freq.hpp \
						utility.cpp  utility.hpp mtrand.cpp \
						figure.cpp figure.hpp \
						seg-site.cpp seg-site.hpp 

hybrid_Lambda_SOURCES = $(shared_SOURCES) main.cpp  
hybrid_Lambda_dbg_SOURCES = $(shared_SOURCES) main.cpp  
hybrid_Lambda_prof_SOURCES = $(shared_SOURCES) main.cpp  

unit_tests_SOURCES = figure.cpp utility.cpp net.cpp node.cpp tests/test_figure.cpp tests/test_runner.cpp

unit_tests_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x 
unit_tests_LDADD= -lcppunit -ldl

include_HEADERS = all_gene_topo.hpp seg-site.hpp sim_gt.hpp freq.hpp utility.hpp mtrand.h param.hpp node.hpp net.hpp


