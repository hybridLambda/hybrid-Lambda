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
#include"net.hpp"
#include<iostream> // clog

enum FIGURE_OPTION { PLOT_DEFAULT, BRANCH, LABEL};
enum FIGURE_PROGRAM { NO_METHOD, LATEX, DOT };

class Figure{
    public:
        //Figure();
        Figure ( int argc, char *argv[] );
        void plot( string net_str );
        
    private:
        FIGURE_PROGRAM method;
        FIGURE_OPTION option;
        void init();
        void check_option();
        void check_method();
        void finalize();
        string figure_file_prefix;
        string figure_file_suffix;
        string figure_file_name;
        Net obj_net;
        valarray <int>  det_x_node ( );
        void plot_in_latex_core();
        void plot_in_latex_file();
        void plot_in_dot( );
        void execute_dot(string method, string suffix);
};
