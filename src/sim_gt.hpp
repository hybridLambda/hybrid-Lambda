/*
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010 -- 2014 Sha (Joe) Zhu
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

/*! \file sim_gt.hpp
 * \brief Header file for sim_gt.cpp */



//#include"utility.hpp"
//#include<math.h>
#include<boost/math/special_functions/gamma.hpp>
#include<boost/math/special_functions/binomial.hpp>
#include<stdio.h>
#include<ctime>
#include"net.hpp"

#ifndef GLOBAL_sim
#define GLOBAL_sim
	


namespace sim{
	class param{
		public:
			param(int argc, char *argv[]);
			
			int num_sim_gt;
			double mutation_rate;
			//double pop_size;
			//double mm;
			//bool pop_size_string_bool;
			bool mm_bool;
			bool pop_bool;
			string sp_string_coal_unit;
			string sp_string_pop_size;
			string para_string;
			vector < int > sample_size;
            string net_str;
			
			//,action_board my_action
			bool num_gener_bool;
			bool sp_coal_unit_bool;
		private:
			param();
	};
}


class action_board{
	public:

	bool sim_mut_unit_bool;
	bool sim_num_gener_bool;
	bool sim_num_mut_bool;
	bool mono_bool;
	bool Si_num_bool;
	//bool total_brchlen_bool;
	string gene_tree_file;
	
	action_board();
	action_board(int argc, char *argv[]);
};


/*! \brief One simulated gene tree from a network under Kingman or multi merger coalescent process*/
class sim_one_gt{
	private:
	void compute_monophyly_vec(Net my_gt_coal_unit,vector < int > sample_size);
	void build_gt_string_mut_unit(double mutation_rate);
	void Si_num_out_table(Net mt_tree,int total_mut);

	
	
	public:
	string gt_string_coal_unit;
	string gt_string_mut_num;
	string gt_string_mut_unit;
	string gt_string_gener_num;
	
	
	vector <string> tax_name;
	vector <double> monophyly;
	double total_brchlen;
	
	sim_one_gt(){
		string gt_string_coal_unit;
		string gt_string_mut_num;
		string gt_string_mut_unit;
		string gt_string_gener_num;
		vector <double> monophyly;
		vector <string> tax_name;
	}
	
	
	//sim_one_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string,vector < int > sample_size, double mutation_rate,bool sim_mut_unit_bool, bool sim_num_gener_bool,bool sim_num_mut_bool,bool mono_bool);
	//sim_one_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string,vector < int > sample_size, double mutation_rate,action_board my_action);
	sim_one_gt(sim::param sim_param,action_board my_action);
	
	//sim_one_gt(string Net_string, vector < int > sample_size, bool multi_merge_bool,double multi_merge_para);
	//void clear(){
		//gt_string_coal_unit.clear();
		//~monophyly;
		//tax_name.clear();
		//}
		
	void clear(){
		gt_string_coal_unit.clear();
		gt_string_mut_num.clear();
		gt_string_mut_unit.clear();
		gt_string_gener_num.clear();
		monophyly.clear();
		tax_name.clear();
		//total_brchlen;
		}
		
};


/*! \brief Collection of simulated gene trees from a network under Kingman or multi merger coalescent process*/
class sim_n_gt{
	public:
		vector <string> gt_string_coal_unit_s;
		vector <string> gt_string_mut_num_s;
		//vector <string> gt_string_gener_num_s;
		//vector <string> gt_string_mut_unit_s;
		
		vector <double> monophyly;
		vector <string> tax_name;
		vector <double> total_brchlen;
		
		sim_n_gt(){
			vector <string> gt_string_coal_unit_s;
			vector <string> gt_string_mut_num_s;
			//vector <string> gt_string_gener_num_s;
			//vector <string> gt_string_mut_unit_s;
			vector <double> monophyly;
			vector <string> tax_name;
			vector <double> total_brchlen;
	
		}
		
		//sim_n_gt(string Net_string, int num_sim_gt,vector <int> sample_size,bool multi_merge_bool,double multi_merge_para);
		//sim_n_gt(string Net_string,int num_sim_gt, string para_string,vector < int > sample_size,double mutation_rate);
		//sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,bool sim_mut_unit_bool, bool sim_num_gener_bool,bool sim_num_mut_bool,bool mono_bool);
		//sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,action_board my_action);
		sim_n_gt(sim::param sim_param,action_board my_action);
		void clear(){
			gt_string_coal_unit_s.clear();
			gt_string_mut_num_s.clear();
			//gt_string_gener_num_s.clear();
			//gt_string_mut_unit_s.clear();
			monophyly.clear();
			tax_name.clear();
			total_brchlen.clear();
		}
};



double Beta(double x,double y);
double unifRand();
int poisson_rand_var(double lambda);
void create_new_site_data(string gt_string_mut_num,int site_i);
string write_sp_string_in_coal_unit(string sp_num_gener_string,string pop_size_string);
string rewrite_pop_string_by_para_string(string para_string,string pop_size_string);
void outtable_header(int total_lineage);
void create_site_data_dir(vector <string> mt_tree_str_s);
void append_seed_to_log_file(unsigned int seed);
valarray <double> build_nc_X(vector < vector <double> > lambda_bk_mat, double num_lineage);
vector < vector <double> > build_lambda_bk_mat(double para,double num_lineage);
int update_nc(valarray <double> nc_X);
double update_coal_para(vector < vector <double> > lambda_bk_mat, double num_lineage);
double kingman_bl(double num_lineage);

#endif
