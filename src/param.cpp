// parameters
#include<param.hpp>


hybridLambda::param::param(){
	seed=(unsigned)(time(0));
	simulation_bool=false;
	help=false;
	freq_bool=false;
	print_tree=false;
	plot_bool;
	log_bool=false;
	log_NAME="LOG";
	seg_bool=false;

};

hybridLambda::param::param(int argc, char *argv[]){
	seed=(unsigned)(time(0));
	simulation_bool=false;
	help=false;
	freq_bool=false;
	print_tree=false;
	plot_bool=false;
	log_bool=false;
	log_NAME="LOG";
	seg_bool=false;	int argc_i=1;
	while (argc_i < argc){
		std::string argv_i(argv[argc_i]);
		
		if (argv_i=="-h" || argv_i=="-help"){
			help=true;
			print_help();
		}
		
		if (argv_i=="-sp_coal_unit" || argv_i=="-sp_num_gener" || argv_i=="-spcu" || argv_i=="-spng"){
			simulation_bool=true;
		}

		if (argv_i=="-seed"){
			read_input_to_param<size_t>(argv[argc_i+1],seed);
			argc_i++;
		}
			
		if (argv_i=="-gt"){
			//reproduce_GENE_trees=false;
			gt_tree_str_s=read_input_lines(argv[argc_i+1]);
			argc_i++;

		}

		if (argv_i=="-mt"){/*! read number of mutations site and simulate segregating sites*/
			mt_tree_str_s=read_input_lines(argv[argc_i+1]);
			argc_i++;

		}

		if (argv_i=="-freq"|| argv_i=="-f" || argv_i=="-freq_file"|| argv_i=="-fF"){
			freq_bool=true;
		}
		
		if (argv_i=="-plot" || argv_i=="-plot_file" || argv_i=="-plotF" || argv_i=="-dot" || argv_i=="-dot_file" || argv_i=="-dotF"){
			plot_bool=true;
		}

		if (argv_i=="-seg"){
			seg_bool=true;
			//sim_num_mut_bool=true;
		}
				
		if (argv_i=="-print"){
			print_tree=true;
		}
			argc_i++;
	
	}
		
	//srand(seed);	// initialize gnu seed
	MTRand_closed mt;
	mt.seed(seed);		// initialize mt seed
}





/*! \brief hybrid-Lambda help file*/
void hybridLambda::print_help(){
	cout<<endl;
	cout<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<"*			hybrid-Lambda beta 0.2			*"<<endl;
	cout<<"*			  Author: Joe ZHU			*"<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<endl<<endl;
	hybridLambda::print_option();
	hybridLambda::print_example();
	exit(1);
}

void hybridLambda::print_example(){
		cout<<"Examples:"<<endl;
	cout<<""<<endl;	
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 3 -seed 2 -gF example1"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp_nt1_para -gF example2 -num 2 -mu 0.00003 -sim mut unit -sim num mut"<<endl;
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num N -pop 25000 -sim num gener"<<endl;
	cout<<"hybrid-Lambda -spng '(A:50000,B:50000)r;' -pop '(A:50000,B:50000)r:40000;'"<<endl;
	cout<<"hybrid-Lambda -spcu '((((A:1.1,B:1.1):2.1,a:2.2):1.1,13D:.2):.3,4:.3);' -S 2 4 3 6 5"<<endl;
	cout<<"hybrid-Lambda -spcu '(A:1,B:1)r;' -mm '(A:1.9,B:.2)r:2;' -S 3 4"<<endl;
	cout<<"hybrid-Lambda -spcu trees/7_tax_sp_nt1_para -dot -branch"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp1 -num 1000 -gF GENE_TREE_FILE -f"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp1 -num 1000 -gF GENE_TREE_FILE -fF FRENQUENCY_FILE"<<endl;	
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 1000 -gF GENE -fF OUTPUT"<<endl;	
	cout<<"hybrid-Lambda -gt GENE_coal_unit -f "<<endl;	
	cout<<"hybrid-Lambda -spcu '(A:5,B:5)r;'-mono -num 100 -mm .1 -S 4 4"<<endl;
	cout<<endl;
}

void hybridLambda::print_option(){
	cout<<setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<endl;
	cout<<setw(20)<<"-spcu INPUT"<<"  --  "<<"Input the species network/tree string through command line or a file."<<endl;
	cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in coalescent unit."<<endl;
	cout<<setw(20)<<"-spng INPUT"<<"  --  "<<"Input the species network/tree string through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in number of generation."<<endl;
	cout<<setw(20)<<"-pop INPUT"<<"  --  "<<"Population sizes are defined by a single numerical constant, "<<endl;
	cout<<setw(26)<<" "<<"or a string which specifies the population size on each branch. "<<endl;
	cout<<setw(26)<<" "<<"The string can be input through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"By default, population size 10,000 is used."<<endl;
	cout<<setw(20)<<"-mm INPUT"<<"  --  "<<"Multiple merger parameters are defined by a single numerical constant, "<<endl;
	cout<<setw(26)<<" "<<"or a string which speifies the parameter on each branch. "<<endl;
	cout<<setw(26)<<" "<<"The string can be input through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"By default, Kingman coalescent is used."<<endl;
	cout<<setw(20)<<"-S n1 n2 ..."<<"  --  "<<"Specify the number of samples for each taxon."<<endl;
	cout<<setw(20)<<"-num N"<<"  --  "<<"The number of gene trees will be simulated."<<endl;
	cout<<setw(20)<<"-seed SEED"<<"  --  "<<"User define random SEED"<<endl;
	cout<<setw(20)<<"-mu MU"<<"  --  "<<"User defined constant mutation rate MU. By default mutation rate 0.00005 is used."<<endl;
	cout<<setw(20)<<"-gF FILE [option]"<<"  --  "<<"Specify the filename for simulated gene trees. \"GENE_TREE\" by default"<<endl;
	//cout<<"     By default, gene tree branch lengths are in coalescent unit "<<endl;
	cout<<setw(20)<<"-sim_mut_unit"<<"  --  "<<"Convert the simulated gene tree branch lengths to mutation unit."<<endl;
	cout<<setw(20)<<"-sim_num_gener"<<"  --  "<<"Convert the simulated gene tree branch lengths to number of generations."<<endl;
	cout<<setw(20)<<"-sim_num_mut"<<"  --  "<<"Simulate numbers of mutations on each branch of simulated gene trees."<<endl;
	cout<<setw(20)<<"-sim_Si_num"<<"  --  "<<"Generate the file out table, which includes the number of segregating"<<endl;
	cout<<setw(26)<<" "<<"sites and the total branch length of the gene tree in coalescent unit."<<endl;
	cout<<setw(20)<<"-f"<<"  --  "<<"Generate a topology frequency table of a set of input trees or simulated gene trees."<<endl;
	cout<<setw(26)<<" "<<"Frequency table is saved in file freq out by default."<<endl;
	cout<<setw(20)<<"-fF FILE"<<"  --  "<<"The topology frequency table will be saved in the FILE."<<endl;
	cout<<setw(20)<<"-gt FILE"<<"  --  "<<"Specify the FILE of trees to analyse tree topology frequencies."<<endl;
	cout<<setw(20)<<"-mono"<<"  --  "<<"Generate a frequency table of monophyletic, paraphyletic and polyphyletic trees. "<<endl;
	cout<<setw(20)<<"-plot/-dot [option]"<<"  --  "<<"Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree)."<<endl;
	cout<<setw(20)<<"      -branch"<<"  --  "<<"Branch lengths will be labelled in the figure."<<endl;
	cout<<setw(20)<<"-plotF/-dotF FILE"<<"  --  "<<"Generated figure will be saved in FILE."<<endl;			
	cout<<endl;	
}
