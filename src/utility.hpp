/* 
 * 
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012, 2013 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-Lambda 
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


/*! \file utility.hpp
 * \brief Header file for network.cpp */

#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<valarray>
#include<cassert>
#include<stdexcept>

using namespace std;

#ifndef GLOBAL_H
#define GLOBAL_H	



/*! 
 */
 

template<class T>
void read_input_to_param(char inchar[], T &input)
{
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
            throw std::invalid_argument("Invalid argument type. ");
	}	
}


bool start_of_tax_name(string in_str,size_t i);


//vector <string> all_n_tax_gene_tree(unsigned int tax_num);
//void plot_in_latex(const char* file_name, string net_str, int plot_option);
//void plot_in_latex_file(const char* file_name, string net_str, int plot_option);
//void plot_in_dot(const char* file_name, string net_str, int plot_option);
//double factorial(double a);	
//double n_permu_a(double n, double a);
//double n_choose_k(double n, double k);
//int factorial_int(int a);	
//int n_permu_a_int(int n, int a);
//int n_choose_k_int(int n, int k);
void appending_debug_file(string debug_file_input);
//void appending_log_file(string log_file_input);/*!< \todo change string to char* type ??*/ 
void appending_log_file(std::string log_file_NAME,std::string log_file_input /*! Information added*/);

string remove_interior_label(string in_str);
string remove_brchlen(string in_str);
string tree_topo(string in_str);


//valarray <int>  det_x_node (Net net_dummy);
//void print_help();
string rm_and_hash_sign(string in_str);
string rm_and_sign(string in_str);
string rm_hash_sign(string in_str);


void check_and_remove(const char* file_name);

int my_exit();


size_t Parenthesis_balance_index_backwards(string in_str,size_t i);
size_t Parenthesis_balance_index_forwards(string in_str,size_t i);
void checking_Parenthesis(string in_str);

string extract_label(string in_str, size_t i);
size_t end_of_label_or_bl(string in_str, size_t i);

string write_para_into_tree(string sp_string, double para);


string extract_label(string in_str, size_t i);
size_t hybrid_hash_index(string in_str);
string extract_hybrid_label(string in_str);
string extract_hybrid_para_str(string in_str);
double extract_hybrid_para(string in_str);

string read_input_line(char inchar[]);
vector <string> read_input_lines(char inchar[]);
string read_input_para(char inchar[],string in_str);
bool is_num(char inchar[]);




/*! \brief Compute factorial of a \return double a! */
template<class T>
T factorial (T a){
	if (a > 1){
		return (a * factorial (a-1));}
	else{
		return (1);}
}


/*! \brief Compute a permutations of n \return double */
template<class T>
T n_permu_a (T n, T a){
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
template<class T>
T n_choose_k(T n, T k){
	if (k<(n/2)){
		return (n_choose_k(n,n-k));}
	else{
		return (n_permu_a(n,k)/factorial(k));}
}




#endif
