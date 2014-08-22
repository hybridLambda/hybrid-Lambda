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


#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <stdio.h>
//#include <ctime>
#include "net.hpp"

#ifndef GLOBAL_sim
#define GLOBAL_sim
	
    
class SimulationParameters{
    friend class HybridLambda;
    friend class sim_one_gt;
    
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
    
    
    bool num_gener_bool;
    bool sp_coal_unit_bool;
    SimulationParameters();
    ~SimulationParameters(){};
    public:
    string net_str;
};


//namespace sim{
	//class param{
		//public:
			//param(int argc, char *argv[]);
			
			//double mutation_rate;
			////double pop_size;
			////double mm;
			////bool pop_size_string_bool;
			//bool mm_bool;
			//bool pop_bool;
			//string sp_string_coal_unit;
			//string sp_string_pop_size;
			//string para_string;
			//vector < int > sample_size;
            //string net_str;
			
			////,action_board my_action
			//bool num_gener_bool;
			//bool sp_coal_unit_bool;
		//private:
			//param();
	//};
//}


class action_board {
    friend class HybridLambda;
    friend class sim_one_gt;
    public:
	bool mono()          const { return mono_bool; }  // \todo, make this private

    private:
	bool sim_mut_unit()  const { return sim_mut_unit_bool;  }
	bool sim_num_gener() const { return sim_num_gener_bool; }
	bool sim_num_mut()   const { return sim_num_mut_bool;   }
	bool Si_num()        const { return Si_num_bool; }

    
    void set_sim_mut_unit()  { this->sim_mut_unit_bool  = true; }
    void set_sim_num_gener() { this->sim_num_gener_bool = true; }
    void set_sim_mut_mut()   { this->sim_num_mut_bool   = true; }	
	bool set_Si_num() {this->Si_num_bool = true; this->sim_num_mut_bool=true; }
    bool set_mono() { this->mono_bool = true; } 
    
	bool sim_mut_unit_bool;
	bool sim_num_gener_bool;
	bool sim_num_mut_bool;
	bool mono_bool;
	bool Si_num_bool;
    
	action_board();
    ~action_board(){};
};


/*! \brief One simulated gene tree from a network under Kingman or multi merger coalescent process*/
class sim_one_gt{
    friend class HybridLambda;
	//private:
        void compute_monophyly_vec(Net my_gt_coal_unit,vector < int > sample_size);
        void build_gt_string_mut_unit(double mutation_rate);
        void Si_num_out_table(Net mt_tree,int total_mut);
        action_board* simulation_jobs_;
        SimulationParameters* parameters_;

	//public:
        string gt_string_coal_unit;
        string gt_string_mut_num;
        string gt_string_mut_unit;
        string gt_string_gener_num;
        
        
        vector <string> tax_name;
        vector <double> monophyly;
        double total_brchlen;
        
    sim_one_gt( SimulationParameters* sim_param, action_board *simulation_jobs);    
    ~sim_one_gt(){};
		
	//void clear(){
		//gt_string_coal_unit.clear();
		//gt_string_mut_num.clear();
		//gt_string_mut_unit.clear();
		//gt_string_gener_num.clear();
		//monophyly.clear();
		//tax_name.clear();
		//}
		
};


double Beta(double x, double y);
double unifRand();
int poisson_rand_var(double lambda);
void create_new_site_data(string gt_string_mut_num,int site_i);
string write_sp_string_in_coal_unit(string sp_num_gener_string,string pop_size_string);
string rewrite_pop_string_by_para_string(string para_string,string pop_size_string);
void outtable_header(int total_lineage);
void create_site_data_dir(vector <string> mt_tree_str_s);
void append_seed_to_log_file(unsigned int seed);
valarray <double> build_nc_X(vector < vector <double> > lambda_bk_mat, double num_lineage);
vector < vector <double> > build_lambda_bk_mat(double para, double num_lineage);
int update_nc(valarray <double> nc_X);
double update_coal_para(vector < vector <double> > lambda_bk_mat, double num_lineage);
double kingman_bl(double num_lineage);

#endif
