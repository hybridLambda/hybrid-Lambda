//parameters header file
#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<valarray>
#include<stdexcept>
#include"mtrand.h"



#ifndef HYBRDRIDLAMBDA_PARAM_INCLUDED
#define HYBRDRIDLAMBDA_PARAM_INCLUDED
using namespace std;



namespace sim{
	class param{
		public:
		param();
		param(int argc, char *argv[]);
		
		string sp_string_coal_unit;
		string sp_string_pop_size;
		string para_string;
		vector < int > sample_size;
		double mutation_rate;
		//,action_board my_action


	};
}

namespace freq{
	class param{
		public:
			bool gene_freq_bool;
			string gene_tree_file;
			bool reproduce_GENE_trees;
	string freq_file_name;

		
		param();
		param(int argc, char *argv[]);
		private:
	};
}




namespace hybridLambda{
	class param{
		public:
				
			size_t seed;			
			bool help;
			bool log_bool;
			string log_NAME;
			
			param();
			param(int argc, char *argv[]);
		
		private:
		
	};
	
	
	void print_help();
	void print_example();
	void print_option();	
}

template<class T>
void read_input_to_param(char inchar[], T &input)
{
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
            throw std::invalid_argument("Invalid argument type. ");
	}	
}

#endif //HYBRDRIDLAMBDA_PARAM_INCLUDED
