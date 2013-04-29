/*
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012, 2013 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-Lambda.
 * 
 * hybrid-Lambda is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/*! \file main.cpp
 *  \brief Main function of hybrid-Lambda */


#include"sim_gt.hpp"
#include"freq.hpp"
#include"figure.hpp"
#include"param.hpp"

using namespace hybridLambda;
//using namespace std;

bool debug_bool;


bool reproduce_GENE_trees;
string gene_tree_file;

//bool sites_data_bool;
//string freq_file_name;

//cstr = new char [str.size()+1];

/*! \todo write a class for all the input*/



int main(int argc, char *argv[]){
	remove("seg-sites");
	check_and_remove("debug_file");
	check_and_remove("log_file");
	//check_and_remove("texfigure.tex");
	//check_and_remove("Total_lengths");

	//const char * freq_file_str;
	vector < int > sample_size;
	debug_bool=false;
	bool sim_n_gt_bool=false;
	
	string net_str;
	//ifstream gt_tree_file;
	//ifstream mt_tree_file;
	//ifstream Net_file;
	//ifstream para_net_file;
	//ifstream pop_size_file;
	
	vector <string> gt_tree_str_s;
	vector <string> mt_tree_str_s;
	//string gt_tree_str;
	//string mt_tree_str;

	int num_sim_gt;
	
	
	string para_string;
	bool para_string_bool=false;
	bool multi_merge_bool=false;
	bool mm_bool=false;
	double multi_merge_para;

	string pop_size_string;
	double pop_size;//=10000;
	bool pop_bool=false;
	
	bool pop_size_string_bool=false;
	bool pop_size_bool=false;
		
	
	unsigned int seed;
	bool seed_bool=false;
	bool print_tree=false;
	bool num_gener_bool=false;

	bool sites_data_bool=false;
	bool sp_coal_unit_bool=false;

	double mutation_rate;//=0.00005;
	bool mutation_rate_bool=false;

	


	bool samples_bool=false;
	
	action_board my_action;
		
	
	
	int total_lineage=0;
	for (int argc_i=0;argc_i<argc;argc_i++){
		string argv_i(argv[argc_i]);
		if (argv_i=="-gt"){
			reproduce_GENE_trees=false;
			gt_tree_str_s=read_input_lines(argv[argc_i+1]);
		}

		if (argv_i=="-mt"){/*! read number of mutations site and simulate segregating sites*/
			gt_tree_str_s=read_input_lines(argv[argc_i+1]);

		}

		if (argv_i=="-sp_coal_unit" || argv_i=="-sp_num_gener" || argv_i=="-spcu" || argv_i=="-spng"){
			if (argv_i=="-sp_num_gener" || argv_i=="-spng"){
				num_gener_bool=true;
			}
			if (argv_i=="-sp_coal_unit" || argv_i=="-spcu"){
				sp_coal_unit_bool=true;
			}
			if (sp_coal_unit_bool && num_gener_bool){
				cout<<"Species tree branch length should only be in Coalescent unit or number of generations. Choose either -sp_num_gener or -sp_coal_unit"<<endl;
				return my_exit();
			}
			net_str=read_input_line(argv[argc_i+1]);
		}


		
		if (argv_i=="-print"){
			print_tree=true;
		}
		
		if (argv_i=="-sim_mut_unit"){
			my_action.sim_mut_unit_bool=true;
		}
		if (argv_i=="-sim_num_gener"){
			my_action.sim_num_gener_bool=true;
		}
		if (argv_i=="-sim_num_mut"){
			my_action.sim_num_mut_bool=true;
		}
		if (argv_i=="-sim_Si_num"){
			my_action.Si_num_bool=true;
			my_action.sim_num_mut_bool=true;
			check_and_remove("out_table");
		}
		
		
		
		if (argv_i=="-seg"){
			sites_data_bool=true;
			my_action.sim_num_mut_bool=true;
		}

		if (argv_i=="-debug"){
			debug_bool=true;
		}
		
		

		
		
		if (argv_i=="-num"){
			sim_n_gt_bool=true;
			string s(argv[argc_i+1]);
			istringstream num_sim_gt_str(s);
			num_sim_gt_str>>num_sim_gt;
		}
		
		if (argv_i=="-seed"){
			string s(argv[argc_i+1]);
			istringstream seed_str(s);
			seed_str>>seed;
			seed_bool=true;
		}

		
		
		

		

		if (argv_i=="-GENE" || argv_i=="-gF"){
			gene_tree_file=argv[argc_i+1];
		}
		
		if (argv_i=="-multi"){
			multi_merge_bool=true;
			string s(argv[argc_i+1]);
			istringstream multi_merge_para_str(s);
			multi_merge_para_str>>multi_merge_para;
			para_string=write_para_into_tree(net_str, multi_merge_para);
		}

		if (argv_i=="-para"){
			para_string_bool=true;
			//para_net_file.open(argv[argc_i+1]);
			//getline (para_net_file,para_string);
			//para_net_file.close();
			para_string=read_input_line(argv[argc_i+1]);
		}
		
		if (argv_i=="-mm"){
			mm_bool=true;
			//string s(argv[argc_i+1]);
			//if (isdigit(s[0])){
				//istringstream multi_merge_para_str(s);
				//double multi_merge_para;
				//multi_merge_para_str>>multi_merge_para;
				//para_string=write_para_into_tree(net_str, multi_merge_para);
			//}
			//else{
				//if (s[0]=='('){
					//para_string=read_input_line(s.c_str());
				//}
				//else{
					//cout<<"Error: check multi merger parameter."<<endl;	
				//}
			//}
			para_string=read_input_para(argv[argc_i+1],net_str);
		}
		
		if (argv_i=="-pop"){
			pop_bool=true;
			pop_size_string=read_input_para(argv[argc_i+1],net_str);
			//cout<<pop_size_string<<endl;
		}
		
		
		if (argv_i=="-pop_size"){
			pop_size_bool=true;
			string s(argv[argc_i+1]);
			istringstream pop_size_str(s);
			pop_size_str>>pop_size;
			pop_size_string=write_para_into_tree(net_str, pop_size);
		}

		if (argv_i=="-pop_string"){
			pop_size_string_bool=true;
			//pop_size_file.open(argv[argc_i+1]);
			//getline (pop_size_file,pop_size_string);
			//pop_size_file.close();
			pop_size_string=read_input_line(argv[argc_i+1]);

		}

		if (argv_i=="-mu"){
			mutation_rate_bool=true;
			string s(argv[argc_i+1]);
			istringstream mut_rate_str(s);
			mut_rate_str>>mutation_rate;
		}
		

		if (argv_i=="-samples" || argv_i=="-S" ){
			samples_bool=true;
			for (int argc_j=argc_i+1;argc_j<argc;argc_j++){
				string s(argv[argc_j]);
				if (!isdigit(s[0])){
					break;
				}
				istringstream pop_num_str(s);
				int pop_num;
				pop_num_str>>pop_num;
				sample_size.push_back(pop_num);
				total_lineage=total_lineage+pop_num;
				
			}
		}
		if (argv_i=="-mono"){
			my_action.mono_bool=true;
		}
	}
	
	if (argc==1 ){
		//scrm_help help;
		//help.print_help();
		hybridLambda::print_help();
	}	//else, proceed

    try {
		time_t start_time = time(0);
	    hybridLambda::param::param hybrid_para(argc, argv);
	    sim::param::param sim_para(argc, argv);
	    figure::param::param figure_para(argc, argv);
	          time_t end_time = time(0);
      
      std::cout << "Simulation took about " << end_time - start_time 
                << " second(s)" << std::endl;
                
      if (hybrid_para.log_bool){          
				std::ofstream log_file;
				log_file.open (hybrid_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
				log_file << "Simulation took about " << end_time - start_time << " second(s) \n";
				log_file.close();
			}
		int sys=system("cat log_file");		

    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}


