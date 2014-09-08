//useless.cpp


//AC::AC(){
	//pAC=0.0;
	//}
class CAC{
	public:
	vector < unsigned int > alive_lineages; 
	vector < class AC* > child_AC;
	
		unsigned int branch_index;
	double pCAC;
		CAC (){
		//pAC=0.0;
		pCAC=0.0;
	}
};
 
class AC{
	public:
	vector < unsigned int > alive_lineages; 
	vector < unsigned int > tip_starting_alive_lineages;
// 	vector < unsigned int > starting_alive_lineages; should be replaced by
	//vector < AC *> starting_AC;

	unsigned int branch_index;
	//vector < int > num_of_lineages_coal_child_AC;
	//vector < AC* > child_AC;
	vector < class CAC * > child_CAC;
	
	//vector < vector < AC* > > child_AC; /*! \todo !!! change to this, instead of vector < AC* > child_AC, some of AC may have multiple ways to be formed */
	double pAC;

	AC (){
		pAC=0.0;
		//pCAC=0.0;
	}
};



/*! Check and remove files*/
void check_and_remove(const char* file_name){
	ifstream my_file(file_name);
	if (my_file.good())
	{
	  remove(file_name);
	}
}



/*! \brief Compute factorial of a \return double a! */
double factorial (double a){
	if (a > 1){
		return (a * factorial (a-1));}
	else{
		return (1);}
}


/*! \brief Compute a permutations of n \return double */
double n_permu_a (double n, double a){
	if (a>1){
		return (n*n_permu_a(n-1,a-1));
	}
	else{
		if (a==1){
			return (n);
		}
		else{
			return (1);
		}
	}
}

/*! \brief Compute n choose k \return double */
double n_choose_k(double n, double k){
	if (k<(n/2)){
		return (n_choose_k(n,n-k));}
	else{
		return (n_permu_a(n,k)/factorial(k));}
}

/*! \brief Compute factorial of a \return int a! */
int factorial_int (int a){
	if (a > 1){
		return (a * factorial_int (a-1));}
	else{
		return (1);}
}

/*! \brief Compute a permutations of n \return int */
int n_permu_a_int (int n, int a){
	if (a>1){
		return (n*n_permu_a_int(n-1,a-1));}
	else{
		if (a==1){
			return (n);}
		else{
			return (1);}
	}
}

/*! \brief Compute n choose k \return int */
int n_choose_k_int(int n, int k){
	if (k<(n/2)){
		return (n_choose_k_int(n,n-k));}
	else{
		return (n_permu_a_int(n,k)/factorial_int(k));}
}


///*! \brief hybrid-Lambda help file*/
//void print_help(){
	//cout<<"*****************************************************************"<<endl;
	//cout<<"*			hybrid-Lambda beta 0.1			*"<<endl;
	//cout<<"*			  Author: Joe ZHU			*"<<endl;
	//cout<<"*****************************************************************"<<endl;
	//cout<<endl<<endl;
	//cout<<"-h or -help -- Help menu"<<endl;
	//cout<<"-num NUMBER -- Define the number of the gene tree will be simulated"<<endl;
	//cout<<"By default, hybrid-Lambda simulates gene trees under Kingmman Coalescent process, with using "<<endl;
	//cout<<"-mm -- it allows to have multi merger process"<<endl;
	//cout<<"-mm ALPHA -- 1 < ALPHA < 2 "<<endl;
	//cout<<"-mm PSI -- 0 < PSI < 1 "<<endl;
	//cout<<"-mm SPECIES_NETWORK_FILE_para -- Define the location of the SPECIES_NETWORK_FILE_para branch length are coalescent parameters in the branch"<<endl;
	//cout<<"-spcu SPECIES_NETWORK_FILE_coal_unit -- Define the location of the SPECIES_NETWORK_FILE branch length in coalescent unit"<<endl;
	//cout<<"-spng SPECIES_NETWORK_FILE_num_gener -- Define the location of the SPECIES_NETWORK_FILE branch length indicate number of generations in the branch"<<endl;
	//cout<<"-pop SPECIES_NETWORK_FILE_pop_size -- Define the location of the SPECIES_NETWORK_FILE branch length indicates population size in the branch"<<endl;
	//cout<<"-pop -- User define constant population size"<<endl;
	//cout<<"-mu -- User define mutation rate"<<endl;
	//cout<<"-gt GENE_TREE_FILE -- Define the location of the GENE_TREE_FILE"<<endl;
	//cout<<"-gF GENE_TREE_FILE -- Define the out name of file that gene trees will be saved, \"GENE_TREE\" by default"<<endl;
	//cout<<"     By default, gene tree branch lengths are in coalescent unit "<<endl;
	//cout<<"     -sim_mut_unit, gene tree branch lengths are in mutation unit "<<endl;
	//cout<<"     -sim_num_gener, gene tree branch lengths are in number of generations "<<endl;
	//cout<<"     -sim_num_mut, gene tree branch lengths are in number of mutations "<<endl;
	//cout<<"     -sim_Si_num, number of segregating sites "<<endl;
	//cout<<"-f -- Count the frequency of the simulated trees, output is saved in \"freq_out\" by default"<<endl;
	//cout<<"-fF FRENQUENCY_FILE -- Define the name of the file that count the frequency of the simulated trees"<<endl;
	//cout<<"-S NUMBER_OF_LINEAGE_ENTERING_1 NUMBER_OF_LINEAGE_ENTERING_1 ... -- Specify number of lineage entering each taxon "<<endl;
	//cout<<"-mono -- Give frequency of topology of monophyly, paraphyly and polyphyly"<<endl;
	//cout<<"-seed -- User define random seed"<<endl;
	//cout<<"-seg -- To produce segregating site data"<<endl;
	
	
	////cout<<"-debug -- To generate a debug file \"debug_file\""<<endl;
			
	//cout<<endl;
	//cout<<"Examples:"<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 "<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -dot "<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -dotF figure.dot "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plot "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plotF figure.tex "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plot -branch"<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plotF -label figure.tex "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000"<<endl;	
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -num 100 -S 100 100"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mm 1.1 -num 100 -S 100 100"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mm 0.1 -num 100 -S 100 100"<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -gF GENE_TREE_FILE"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -f"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -fF FRENQUENCY_FILE"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -gF GENE_TREE_FILE -f"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -gF GENE_TREE_FILE -fF FRENQUENCY_FILE"<<endl;	
	//cout<<"hybrid-Lambda -gt GENE_TREE_FILE -f "<<endl;	
	//cout<<"hybrid-Lambda -gt GENE_TREE_FILE -fF FRENQUENCY_FILE"<<endl;	
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mono -num 100 -mm .1 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mono -num 100 -seed 2 -mm .1 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -para 2_tax_sp1_para2 -mono -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mm 2_tax_sp1_para2 -seg -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mm 2_tax_sp1_para2 -seg -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_size 10000 -para 2_tax_sp1_para2 -seg -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mu 0.0002 -para 2_tax_sp1_para2 -seg -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mu 0.00002  -seg -num 1000 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mu 0.00002 -num 1000 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_size 10000 -mu 0.0002 -mm 1.4 -seg -num 100 -S 4 4"<<endl;
	//cout<<endl;
//}

void appending_log_file(string log_file_input /*! Information added*/){
	ofstream log_file;
	log_file.open ("log_file", ios::out | ios::app | ios::binary); 
	log_file << log_file_input << "\n";
	log_file.close();
}






		string appending_log_str;

		
		if (reproduce_GENE_trees){
			//string gene_tree_file_s="rm "+gene_tree_file+"*";
			string gene_tree_file_coal_unit=gene_tree_file+"_coal_unit";
			string gene_tree_file_mut_unit=gene_tree_file+"_mut_unit";
			string gene_tree_file_num_gener=gene_tree_file+"_num_gener";
			string gene_tree_file_num_mut=gene_tree_file+"_num_mut";
			remove(gene_tree_file_coal_unit.c_str());
			remove(gene_tree_file_mut_unit.c_str());
			remove(gene_tree_file_num_gener.c_str());
			remove(gene_tree_file_num_mut.c_str());
		}
		else{
			if (net_str.size()>0){
				net_str.clear();
				appending_log_file("Illegal flags");
				cout<<"Programs terminated"<<endl;
				//return 0;
				return my_exit();
			}
		}


				
		if (net_str.size()>0){

						
			if (!samples_bool){
				total_lineage=net_dummy.tax_name.size();
				for (int i=0;i<total_lineage;i++){
					sample_size.push_back(1);
				}
			}
			else{//  check the number of lineages and the number of species 
				if (sample_size.size()!=net_dummy.tax_name.size()){
					appending_log_file("Numbers of samples and numbers of species not equal!!!");
					appending_log_file("Simulation terminated");
					return my_exit(); 
				}
			}	
			
			//if ((!multi_merge_bool) && (!para_string_bool)){
			//if (!mm_bool){
			if ((!multi_merge_bool) && (!para_string_bool) && !mm_bool){

				para_string=write_para_into_tree(net_str, 2.0); // If coalescent parameter is ungiven, use Kingman coalescent as default
				appending_log_file("Default Kingman coalescent on all branches");
			}
	
			//if ((!pop_size_bool) && (!pop_size_string_bool)){
			//if (!pop_bool){
			if ((!pop_size_bool) && (!pop_size_string_bool) && (!pop_bool)){
				pop_size_string=write_para_into_tree(net_str, 10000.0); // If the population size is ungiven, use default population size 10000
				if (my_action.sim_num_gener_bool || my_action.sim_mut_unit_bool || my_action.sim_num_mut_bool || my_action.Si_num_bool){
					appending_log_file("Default population size of 10000 on all branches");
					
				}
			}
			
			//appending_log_file(pop_size_string);
			
			//if (!mutation_rate_bool){
				//mutation_rate=0.00005;
				//if (my_action.sim_mut_unit_bool || my_action.sim_num_mut_bool || my_action.Si_num_bool){
					//appending_log_file("Default mutation rate 0.00005");
				//}
			//}
			
			pop_size_string=rewrite_pop_string_by_para_string(para_string,pop_size_string);  // checking if modify pop_size_string is needed,
	
			if (num_gener_bool){
				net_str=write_sp_string_in_coal_unit(net_str,pop_size_string);	// Convert number of generations and population size to coalescent unit
			}
			
			
			Net new_net_dummy(net_str);
			
			if (!new_net_dummy.is_ultrametric){
				appending_log_file("WARNING! NOT ULTRAMETRIC!!!");
				//return my_exit();
			}
		
			append_seed_to_log_file(seed);

			//if (!sim_n_gt_bool){
				//num_sim_gt=1;
			//}
			if (my_action.Si_num_bool){
				outtable_header(total_lineage);
			}
			//cout<<net_str<<endl;
			//cout<<pop_size_string<<endl;
			sim_n_gt simd_gt_tree_str_s(net_str,pop_size_string,para_string,sample_size,mutation_rate,num_sim_gt,my_action);
			gt_tree_str_s=simd_gt_tree_str_s.gt_string_coal_unit_s;
			mt_tree_str_s=simd_gt_tree_str_s.gt_string_mut_num_s;
			
			if (my_action.mono_bool  && sample_size.size()==2){
				cout<<"   A mono     B mono Recip mono     A para     B para  Polyphyly"<<endl;
				for (unsigned int mono_i=0;mono_i<simd_gt_tree_str_s.monophyly.size();mono_i++){
					cout<<setw(9)<<simd_gt_tree_str_s.monophyly[mono_i]<<"  ";
				}
				cout<<endl;
			}
			
			ostringstream num_ostr_stream;
			num_ostr_stream<<num_sim_gt;
			appending_log_str=num_ostr_stream.str() + " trees simulated.";
			appending_log_file(appending_log_str);
			
			if (!sim_n_gt_bool){
				string command="cat "+gene_tree_file+"_coal_unit";
				int sys=system(command.c_str());
			}
			
		}
		
		if (my_action.Si_num_bool){
			appending_log_file("Table of number of segregating site in file: out_table");
		}
		
		if (sites_data_bool){
			create_site_data_dir(mt_tree_str_s);
		}
		
		if (gene_freq_bool){
			compute_gt_frequencies(gt_tree_str_s, freq_file_name);
		}
		


void appending_debug_file(string debug_file_input){
	ofstream debug_file;
	debug_file.open ("debug_file", ios::out | ios::app | ios::binary); 
	debug_file << debug_file_input << "\n";
	debug_file.close();
}



/*! \brief Record the used random seed into the log_file*/
void append_seed_to_log_file(unsigned int seed){
	ostringstream seed_ostr_stream;
	seed_ostr_stream<<seed;
	string appending_log_str="Random Seed  " + seed_ostr_stream.str() + "  used ";
	//appending_log_file(appending_log_str);
}




/*! \brief Add more information to log_file */
void appending_log_file(std::string log_file_NAME,std::string log_file_input /*! Information added*/){
	std::ofstream log_file;
	log_file.open (log_file_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
	log_file << log_file_input << "\n";
	log_file.close();
}


/*! \brief Terminate the program and print out the log file*/
int my_exit(){ // change my_exit() to return 0, and terminate program
	int sys=system("cat log_file");		
	//exit(1);
	return 0;
}

/*! \brief UNUSED AT THE MOMENT. Add new simulated gene trees into the list of gene trees 
 * \todo change the target gene tree file name, */
void appending_sim_gt_file(string gene_tree_file,string sim_gt_input){
	ofstream sim_gt_file;
	sim_gt_file.open (gene_tree_file.c_str(), ios::out | ios::app | ios::binary); 
	sim_gt_file << sim_gt_input << "\n";
	sim_gt_file.close();
}
