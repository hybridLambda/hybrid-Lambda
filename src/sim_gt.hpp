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
#include "net.hpp"

#ifndef GLOBAL_sim
#define GLOBAL_sim
	
class SimulationParameters{
    friend class HybridLambda;
    friend class sim_one_gt;
    
    double mutation_rate;
    bool mm_bool;
    bool pop_bool;
    bool samples_bool;
    bool is_Net;
    string sp_string_coal_unit;
    string sp_string_pop_size;
    string para_string;
    int total_num_lineage;
    
    bool num_gener_bool;
    bool sp_coal_unit_bool;
    SimulationParameters();
    ~SimulationParameters(){};
    void finalize();
    string net_str;
    vector < int > sample_size;
};


class action_board {
    friend class HybridLambda;
    friend class sim_one_gt;

    bool sim_mut_unit()  const { return sim_mut_unit_bool;  }
    bool sim_num_gener() const { return sim_num_gener_bool; }
    bool sim_num_mut()   const { return sim_num_mut_bool;   }
    bool Si_num()        const { return Si_num_bool; }    
    
    void set_sim_mut_unit()  { this->sim_mut_unit_bool  = true; }
    void set_sim_num_gener() { this->sim_num_gener_bool = true; }
    void set_sim_num_mut()   { this->sim_num_mut_bool   = true; }	
    bool set_Si_num() { this->Si_num_bool = true; this->sim_num_mut_bool=true; }
    bool set_mono() { this->mono_bool = true; } 
    bool mono()          const { return mono_bool; }  // \todo, make this private
    
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

    void compute_monophyly_vec(Net &my_gt_coal_unit,vector < int > sample_size);
    void Si_num_out_table(Net &mt_tree);
    action_board* simulation_jobs_;
    SimulationParameters* parameters_;
    string gt_string_coal_unit;
    string gt_string_mut_num;
    string gt_string_mut_unit;
    string gt_string_gener_num;
    vector <string> tax_name;
    vector <double> monophyly;
    double total_brchlen;
    int total_mut;    
    ofstream * Si_table_;
    sim_one_gt( SimulationParameters* sim_param, action_board *simulation_jobs, std::ofstream &Si_table );    
    ~sim_one_gt(){};		
    
    vector < vector <double> > build_lambda_bk_mat(double para, double num_lineage);
    void build_gt_string_mut_unit();
    void sim_mt_tree();


};


double Beta(double x, double y);
double unifRand();
int poisson_rand_var(double lambda);
string write_sp_string_in_coal_unit(string sp_num_gener_string, string pop_size_string);
string rewrite_pop_string_by_para_string(string para_string,string pop_size_string);

valarray <double> build_nc_X(vector < vector <double> > lambda_bk_mat, double num_lineage);
int update_nc(valarray <double> nc_X);
double update_coal_para(vector < vector <double> > lambda_bk_mat, double num_lineage);
double kingman_bl(double num_lineage);

string write_para_into_tree(string sp_string, double para);

string construct_adding_new_Net_str(Net & old_Net);

#endif
