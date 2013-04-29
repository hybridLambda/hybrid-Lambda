// parameters
#include<param.hpp>


hybridLambda::param::param(){
	
};

hybridLambda::param::param(int argc, char *argv[]){
	param();
		//if (!seed_bool){
			//seed=(unsigned)(time(0));
		//}

		//srand(seed);	// initialize gnu seed
		MTRand_closed mt;
		mt.seed(seed);		// initialize mt seed
		
		int argc_i=1;
		while (argc_i < argc){
					if (argv_i=="-freq"|| argv_i=="-f"){
				gene_freq_bool=true;
			}
			if (argv_i=="-freq_file"|| argv_i=="-fF"){
				gene_freq_bool=true;
				freq_file_name=argv[argc_i+1];
			}
			check_and_remove(freq_file_name.c_str());	
			
		}
		
		
		
}


sim::param::param(){
	seed=(unsigned)(time(0));
}

sim::param::param(){
	param();
		
		int argc_i=1;
	while( argc_i < argc ){
		
		std::string argv_i(argv[argc_i]);
		
		if (argv_i=="-nsam"){ // if scrm is not called, use this option read in the number of samples
			//read_input_to_int(argv[argc_i+1],nsam);
			read_input_to_param<int>(argv[argc_i+1],nsam);
			argc_i++;
		}
	}
	
	
}

figure::param::param(int argc, char *argv[]){
	
}

figure::param::param(){
		Net net_dummy(net_str);

	if (print_tree){
		net_dummy.print_all_node();
		appending_log_file("Tree printed");
	}

	plot_option=set_plot_option(plot_label,plot_branch);

	if (plot_bool){
		plot_in_latex_file(tex_fig_name.c_str(), net_dummy,plot_option);	
	}
	
	if (dot_bool){
		plot_in_dot(dot_fig_name.c_str(), net_dummy,plot_option);			
	}
	
	if (print_tree || plot_bool || dot_bool){
		return my_exit();
	}
}
