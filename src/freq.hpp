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

/*! \file freq.hpp
 * \brief Header file for freq.cpp */


#include"utility.hpp"
#include"all_gene_topo.hpp"
#include"sim_gt.hpp"


/*! \brief Class function for counting topological frequencies */
class topo_freq{
	public:
	vector <string> gene_topo; /*!< \brief Different gene topologies */
	vector <int> gene_freq;    /*!< \brief Frequencies for different tree topologies */
	topo_freq(vector <string> gt_strings); /*!< \brief topo_freq class constructor */
};

bool same_topo(string gt_string1,string gt_string2);

void compute_gt_frequencies(vector <string> gt_tree_str_s, string freq_file_name);
