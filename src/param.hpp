//parameters header file
#include"mtrand.h"
#include"utility.hpp"

#ifndef HYBRDRIDLAMBDA_PARAM_INCLUDED
#define HYBRDRIDLAMBDA_PARAM_INCLUDED
using namespace std;

namespace hybridLambda{
	class param{
		public:	
			size_t seed;				
			bool simulation_bool;
			bool help;
			bool freq_bool;
			bool print_tree;
			bool plot_bool;
			bool log_bool;
			string log_NAME;
			bool seg_bool;
			//bool mono_bool;		
			//bool sites_data_bool;
			
			vector <string> gt_tree_str_s;
			vector <string> mt_tree_str_s;
			
			param(int argc, char *argv[]);
		
		private:
			param();
		
	};
	
	void print_help();
	void print_example();
	void print_option();	
}

#endif //HYBRDRIDLAMBDA_PARAM_INCLUDED
