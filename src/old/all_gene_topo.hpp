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

/*! \file all_gene_topo.hpp
 * 	\brief Header file for all_gene_topo.cpp */

//#include<iostream>
//#include<vector>
//#include<string>
#include"utility.hpp"

using namespace std;

vector <string> generate_new_topo_list(vector <string> old_topo_list, vector <string> taxa_name);
vector <string> all_n_tax_gene_tree(vector <string> taxa_name);
vector < string > find_current_taxa_name(string in_str);
string add_new_taxa_at_tip(string in_str,size_t i, size_t tax_i, vector <string> taxa_name,string label);
string add_new_taxa_at_int(string in_str, size_t i, size_t tax_i, vector <string> taxa_name);
