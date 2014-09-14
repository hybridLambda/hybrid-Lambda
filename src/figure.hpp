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
#include <iostream> // clog

enum FIGURE_OPTION { PLOT_DEFAULT, BRANCH, LABEL};
enum FIGURE_PROGRAM { NO_METHOD, LATEX, DOT };

class Figure{
    friend class HybridLambda;
    #ifdef UNITTEST
    friend class TestFigure;
    #endif
    // Members
    string figure_file_prefix;
    FIGURE_PROGRAM method;
    FIGURE_OPTION option;    
    int argc_;
    int argc_i;
    char * const* argv_;
    valarray <int> x_node;
    vector <int> x_node_tmp;
    vector <size_t> x_node_tmp_index;
    ofstream figure_ofstream;
    string figure_file_suffix;
    string figure_file_name;
    Net obj_net;
    
    // Methods    
    Figure ( int argc, char * const* argv );
    ~Figure(){}
    void plot( string net_str );
    void init();
    void initialize_method( FIGURE_PROGRAM program, string suffix);
    void check_option();
    void check_method();
    void finalize();
    void det_x_node ( );
    void x_node_shift();
    void plot_in_latex();
    void plot_in_dot( );
    void plot_core();
    void execute_dot(string method, string suffix);
    void edge_entry(string from, string to, size_t label, double bl, bool tip);

    /*! \brief Remove '#' signs and the gamma parameter from a string \return string */
    string rm_hash_sign( string in_str ){
        while ( int(in_str.find('#'))>0 && in_str.find('#')!=string::npos ) {
            size_t i = end_of_label_or_bl(in_str, in_str.find('#'));
            in_str.erase(in_str.find('#'),i-in_str.find('#')+1);
            //in_str.erase(in_str.find('#'),1);
        }
        return in_str;
    }

    /*! \brief Remove the '&' and '#' signs from a string \return string */
    string rm_and_hash_sign( string in_str ){
        return rm_hash_sign( rm_and_sign ( in_str ) );
    }

    string rm_and_sign( string in_str ){
        while ( int(in_str.find('&')) > 0 && in_str.find('&')!=string::npos) {
            in_str.erase(in_str.find('&'),1);
        }
        return in_str;
    }
};
