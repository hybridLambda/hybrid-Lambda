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
//parameters header file
#include"mtrand.h"
#include"utility.hpp"

#ifndef HYBRDRIDLAMBDA_PARAM_INCLUDED
#define HYBRDRIDLAMBDA_PARAM_INCLUDED
using namespace std;

namespace hybridLambda{
	class param{
		public:	
			size_t seed;				
			bool simulation_bool;
			bool help;
			bool freq_bool;
			bool print_tree;
			bool plot_bool;
			bool log_bool;
			string log_NAME;
			bool seg_bool;
			bool read_GENE_trees;
			bool read_mt_trees;
			string gt_file_name;
			string mt_file_name;
			
			bool mm_bool;
			bool pop_bool;
			bool tmrca_bool;
			bool bl_bool;
			string tmrca_NAME;
			string bl_NAME;
			//string net_str;
			//bool mono_bool;		
			//bool sites_data_bool;
			
			//vector <string> gt_tree_str_s;
			//vector <string> mt_tree_str_s;
			
			param(int argc, char *argv[]);
		
		private:
			param();
			void init();
	};
	
	void print_help();
	void print_example();
	void print_option();	
}

#endif //HYBRDRIDLAMBDA_PARAM_INCLUDED
