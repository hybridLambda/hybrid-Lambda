//parameters header file


#ifndef HYBRDRIDLAMBDA_PARAM_INCLUDED
#define HYBRDRIDLAMBDA_PARAM_INCLUDED



namespace sim{
	class param{
public:
param(int argc, char *argv[]);

string sp_string_coal_unit;
string sp_string_pop_size;
string para_string;
vector < int > sample_size;
double mutation_rate;
//,action_board my_action


	}
}

namespace figure{
	class param{

	}
}



namespace hybridLambda{
	class param{
		public:
		
		log_NAME="scrm.log";

		
		}
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
