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

//net.hpp
#include "node.hpp"
#include <valarray>

/*!\file net.cpp
 *  \brief Core function of converting a Newick (extended Newick) format string into a species tree (network), and simple string manipulation for tree strings
*/

#ifndef NETWORK
#define NETWORK


/*! \brief Network class*/
class Net{
    friend class HybridLambda;
    friend class sim_one_gt;
	private:
		
		string label_interior_node(string in_str);
		void enumerate_internal_branch( Node &node );
		
	    size_t first_coal_rank();
        size_t current_enum_;
        void init(){
            this->current_enum_ = 0;
            this->is_Net = false;
            this->is_ultrametric = true;
            }

        bool start_of_tax_name(string in_str, size_t i);
        size_t Parenthesis_balance_index_backwards( string &in_str, size_t i );
        size_t Parenthesis_balance_index_forwards( string &in_str, size_t i );
        
        void check_Parenthesis(string &in_str);
        void check_labeled( string in_str );
        
        void check_isNet(); /*!< \brief To determin if a Net is network or not. \return is_Net */
		void check_isUltrametric(); /*!< \brief To determin if a Net is ultrametric or not. \return is_ultrametric */

        size_t first_coal_index ();
        //void clear(); 
		void print_all_node();
		bool print_all_node_dout();
        
        string rewrite_internal_node_content( size_t i);       
        void connect_graph();
        void extract_tax_and_tip_names();
        
        void init_descendant();
        void init_node_clade();
        void rewrite_descendant();
        void rewrite_node_clade();
        
	public:	
		void rewrite_node_content();
        string net_str; /*!< \brief species network string \todo this is new!!!*/
		size_t max_rank;
		vector< valarray <int> > descndnt;
		vector< valarray <int> > descndnt2;
		vector<string> tax_name;
		vector<string> tip_name;
		//vector <Node*> NodeContainer_ptr; /*!< \brief pointers to the nodes of Net \todo use this instead of NodeContainer */
		vector <Node> NodeContainer;  /*!< \brief vector of nodes */
		bool is_Net; /*!< \brief true if Net is a network; false if it's a tree */
		bool is_ultrametric; /*!< \brief true if the distances between tips and root are equal; false, otherwise */
		

        ~Net(){};
        Net (){ this->init(); }
		//Net (){
			//string net_str;
			//vector <string> tax_name;
			//vector <string> tip_name;
			//vector <Node> NodeContainer;
			////vector <Node*> NodeContainer_ptr;
		//}
				
		Net(string Net_str);
};

#endif
