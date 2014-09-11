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
#include "mtrand.h"
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

    string net_str;
    vector < int > sample_size;

    SimulationParameters();
    ~SimulationParameters(){};
    string write_sp_string_in_coal_unit( string &sp_num_gener_string, string &pop_size_string );
    string rewrite_pop_string_by_para_string( string para_string,string pop_size_string );

    void finalize();
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

    action_board* simulation_jobs_;
    SimulationParameters* parameters_;
    string gt_string_coal_unit;
    string gt_string_mut_num;
    string gt_string_mut_unit;
    string gt_string_gener_num;
    //vector <string> tax_name;
    vector <double> monophyly;
    double total_brchlen;
    int total_mut;    
    ofstream * Si_table_;

    sim_one_gt( SimulationParameters* sim_param, action_board *simulation_jobs, std::ofstream &Si_table );    
    ~sim_one_gt(){};		

	Net my_gt_coal_unit;
	Net my_gt_num_gener;
    void init();
bool sim_num_gener_bool_;
	vector <size_t> remaining_sp_node;
void initialize_gt_tip_nodes( Net & my_Net );

Node *current_sp_pop_node;

    void compute_monophyly_vec( Net &my_gt_coal_unit,vector < int > sample_size );
    void Si_num_out_table ( Net &mt_tree );
    
    vector < vector <double> > build_lambda_bk_mat(double para, double num_lineage);
    void build_gt_string_mut_unit();
    void build_mt_tree();
    
    double update_coal_para( vector < vector <double> > &lambda_bk_mat, double num_lineage);
    valarray <double> build_nc_X(vector < vector <double> > &lambda_bk_mat, double num_lineage);
    int update_nc(valarray <double> nc_X);

    string extract_hybrid_para_str(string in_str){
        size_t hash_index = hybrid_hash_index(in_str);
        return in_str.substr(hash_index+1);//,in_str.size()-1);
    }
    
    double extract_hybrid_para( string in_str ){
        double para;
        istringstream para_istr(extract_hybrid_para_str(in_str));
        para_istr>>para;
        return para;
    }

    /*! \brief Beta function, requires tgamma function from math.h \return double */
    double Beta(double x,double y){
        double Beta_return;
    //	Beta_return=tgamma(x)*tgamma(y)/tgamma(x+y);
        Beta_return=exp(log(tgamma(x))+log(tgamma(y))-log(tgamma(x+y)));
        return Beta_return;
    }
    
    /*! \fn double unifRand()
     * \brief Simulate random variable between 0 and 1.
     */
    double unifRand(){
        MTRand_closed return_value;
        return return_value();
        //return rand()/(double(RAND_MAX)+1);
        //return rand() / double(RAND_MAX); //generates a psuedo-random float between 0.0 and 0.999...
    } 
    
    double kingman_bl( double num_lineage ){
	    return -log( 1-unifRand() ) / num_lineage / (num_lineage-1) * 2.0;
    }
    
    /*! \fn int poisson_rand_var(double lambda)
     * \brief Simulating Poisson random variable from given lambda
     */
    int poisson_rand_var( double lambda ){
        double L = exp(-lambda);
        int k = 0;
        double p = 1;
        while ( p > L ){
             k++; //k =k + 1;
             p *= unifRand(); //p=p*unifRand();
        }
        //Generate uniform random number u in [0,1] and let pp × u.
        k--;//k=k-1;	
        return k;
    }
    
    /*! \brief Compute factorial of a \return double a! */
    template < class T > T factorial ( T a ){
        if (a > 1) return (a * factorial (a-1));
        else       return (1);
    }
    
    /*! \brief Compute a permutations of n \return double */
    template < class T > T n_permu_a ( T n, T a ){
        if   ( a > 1 ) return (n*n_permu_a(n-1,a-1));
        else if (a==1) return (n);
        else           return (1);
    }
    
    /*! \brief Compute n choose k \return double */
    template < class T > T n_choose_k ( T n, T k ){
        if ( k < ( n/2 ) ) return (n_choose_k(n,n-k));
        else               return (n_permu_a(n,k)/factorial(k));
    }

};

string write_para_into_tree(string sp_string, double para);
string construct_adding_new_Net_str(Net & old_Net);

#endif
