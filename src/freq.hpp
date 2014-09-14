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

#include "net.hpp"

class Frequency{
    friend class HybridLambda;
    // Methods
    Frequency (){};
    Frequency ( int argc, char * const* argv );
    ~Frequency(){}
    void compute_gt_frequencies( vector <string> &gt_tree_str_s );
    string tree_topo( string in_str );
    string remove_brchlen( string in_str );
    bool same_topo( string gt_string1, string gt_string2 );
    void compute_gt_frequencies_core();
    
    // Members
    string freq_out_filename;
    int argc_;
    int argc_i;
    char * const* argv_;
    ofstream freq_out_file;
    vector <string> gene_topo; /*!< \brief Different gene topologies */
    vector <int> gene_freq;    /*!< \brief Frequencies for different tree topologies */        
    vector <string> gt_tree_str_tmp;
};
