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

//parameters header file
#include"mtrand.h"
#include"sim_gt.hpp"

#ifndef HYBRDRIDLAMBDA_PARAM_INCLUDED
#define HYBRDRIDLAMBDA_PARAM_INCLUDED
using namespace std;

#include <iostream>     // std::cout, std::endl
#include <iomanip>      // std::setw
#include <stdlib.h>     /* exit, EXIT_FAILURE */


void print_example();
void print_help();
void print_option();


class HybridLambda{
	//class param{
    public:	
        /*!
         * Constructors and Destructors
         */  
        //param();
        
        HybridLambda(int argc, char *argv[]) : argc_(argc), argv_(argv) { init(); }
        void parse() ;
        //HybridLambda(int argc, char *argv[]);
        void HybridLambda_core(sim::param sim_param,action_board my_action);

        /*!
         * Members
         */              
        size_t seed;				
        string gt_file_name;
        string mt_file_name;
        
        bool fst_bool;
        void extract_tmrca();
        void extract_bl();
        void extract_firstcoal();

        vector <string> gt_tree_str_s;
        vector <string> mt_tree_str_s;

        bool freq_bool;
        bool print_tree;
        bool plot_bool;

        bool seg_bool;
        bool read_GENE_trees;
        bool read_mt_trees;
        bool simulation_bool;
        vector <double> monophyly;

    private:
        action_board* simulation_jobs;
        
        
        void init();
        std::ofstream extract_file;
        string extract_file_name;
        
        // Extract time from can first coalscent from the gene trees
        bool tmrca_bool;
        bool bl_bool;        
        bool firstcoal_bool;
        
       	ofstream sim_gt_file_coal_unit;
        ofstream sim_gt_file_mut_unit;
        ofstream sim_gt_file_num_gener;
        ofstream sim_gt_file_num_mut;

        string prefix;
        
        bool mm_bool;
        bool pop_bool;
        
		vector <string> tax_name;
          int argc_;
  int argc_i;
  char * const* argv_;
  size_t random_seed_;  
};
	

//}


///*! \brief Collection of simulated gene trees from a network under Kingman or multi merger coalescent process*/
//class sim_n_gt{
	//public:
		//vector <string> gt_string_coal_unit_s;
		//vector <string> gt_string_mut_num_s;
		////vector <string> gt_string_gener_num_s;
		////vector <string> gt_string_mut_unit_s;
		
		//vector <double> monophyly;
		//vector <string> tax_name;
		//vector <double> total_brchlen;
		
		//sim_n_gt(){
			//vector <string> gt_string_coal_unit_s;
			//vector <string> gt_string_mut_num_s;
			////vector <string> gt_string_gener_num_s;
			////vector <string> gt_string_mut_unit_s;
			//vector <double> monophyly;
			//vector <string> tax_name;
			//vector <double> total_brchlen;
	
		//}
		
		////sim_n_gt(string Net_string, int num_sim_gt,vector <int> sample_size,bool multi_merge_bool,double multi_merge_para);
		////sim_n_gt(string Net_string,int num_sim_gt, string para_string,vector < int > sample_size,double mutation_rate);
		////sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,bool sim_mut_unit_bool, bool sim_num_gener_bool,bool sim_num_mut_bool,bool mono_bool);
		////sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,action_board my_action);
		//sim_n_gt(sim::param sim_param,action_board my_action);
		//void clear(){
			//gt_string_coal_unit_s.clear();
			//gt_string_mut_num_s.clear();
			////gt_string_gener_num_s.clear();
			////gt_string_mut_unit_s.clear();
			//monophyly.clear();
			//tax_name.clear();
			//total_brchlen.clear();
		//}
//};

#endif //HYBRDRIDLAMBDA_PARAM_INCLUDED
