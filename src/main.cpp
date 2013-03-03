/*
 * hybrid_sim is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_sim.
 * 
 * hybrid_sim is free software: you can redistribute it and/or modify
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
 *  \brief Main function of hybrid_sim */


#include"sim_gt.hpp"
#include"freq.hpp"
//using namespace std;

bool debug_bool;
//bool reproduce_GENE_trees;
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
	bool help=false;
	bool sim_n_gt_bool=false;
	
	string net_str;
	//ifstream gt_tree_file;
	//ifstream mt_tree_file;
	//ifstream Net_file;
	//ifstream para_net_file;
	//ifstream pop_size_file;
	string freq_file_name="freq_out";
	string tex_fig_name="texfigure.tex";
	string dot_fig_name="figure.dot";
	vector <string> gt_tree_str_s;
	vector <string> mt_tree_str_s;
	//string gt_tree_str;
	//string mt_tree_str;

	int num_sim_gt;
	bool gene_freq_bool=false;
	gene_tree_file="GENE_TREE";
	bool reproduce_GENE_trees=true;
	
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

	
	bool plot_bool=false;
	bool dot_bool=false;
	int plot_option;//=0;
	bool plot_label=false;
	bool plot_branch=false;

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

		if (argv_i=="-h" || argv_i=="-help"){
			help=true;
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
		
		
		if (argv_i=="-freq"|| argv_i=="-f"){
			gene_freq_bool=true;
		}
		if (argv_i=="-freq_file"|| argv_i=="-fF"){
			gene_freq_bool=true;
			freq_file_name=argv[argc_i+1];
		}
		check_and_remove(freq_file_name.c_str());
		
		
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

		if (argv_i=="-plot"){
			plot_bool=true;
		}
		if (argv_i=="-plot_file" || argv_i=="-plotF"){
			plot_bool=true;
			tex_fig_name=argv[argc_i+1];
		}
		check_and_remove(tex_fig_name.c_str());
		
		
		if (argv_i=="-label"){
			plot_label=true;
		}

		if (argv_i=="-branch"){
			plot_branch=true;
		}		

		if (argv_i=="-dot"){
			dot_bool=true;
		}
		if (argv_i=="-dot_file" || argv_i=="-dotF"){
			dot_bool=true;
			dot_fig_name=argv[argc_i+1];
		}
		check_and_remove(dot_fig_name.c_str());
		

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
	
	if (help || argc==1){
		print_help();
	}
	else{
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

		if (!seed_bool){
			seed=(unsigned)(time(0));
		}

		//srand(seed);	// initialize gnu seed
		MTRand_closed mt;
		mt.seed(seed);		// initialize mt seed
				
		if (net_str.size()>0){
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
			
			if (!mutation_rate_bool){
				mutation_rate=0.00005;
				if (my_action.sim_mut_unit_bool || my_action.sim_num_mut_bool || my_action.Si_num_bool){
					appending_log_file("Default mutation rate 0.00005");
				}
			}
			
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

			if (!sim_n_gt_bool){
				num_sim_gt=1;
			}
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
		
		int sys=system("cat log_file");		
	}
	return 0;
}



