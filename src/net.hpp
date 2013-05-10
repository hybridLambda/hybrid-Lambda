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
//net.hpp
//#include"utility.hpp"
#include"node.hpp"
//#include<cassert>

/*!\file net.cpp
 *  \brief Core function of converting a Newick (extended Newick) format string into a species tree (network), and simple string manipulation for tree strings
*/

#ifndef NETWORK
#define NETWORK

/*! \brief Collection of Networks in a vector container*/
class Network_s{
	public:
	vector <class Net *> Net_vec;
	void clear(){Net_vec.clear();};
	
	Network_s(){
		vector <class Net *> Net_vec;
	}
	
};

/*! \brief Network class*/
class Net{	
	private:
		string checking_labeled(string in_str);
		string label_interior_node(string in_str);
		int enumerate_internal_branch(Node *current,int e_num_old);
		bool is_net_func(); /*!< \brief To determin if a Net is network or not. \return is_net */
		bool is_ultrametric_func(); /*!< \brief To determin if a Net is ultrametric or not. \return is_ultrametric */
	
	public:
	
		string net_str; /*!< \brief species network string \todo this is new!!!*/
	//	class Network_s SubNetworkS; /*!< \brief sub species networks \todo this is new!!!*/
		int max_rank;
		vector< valarray <int> > descndnt;
		vector< valarray <int> > descndnt2;
		vector<string> tax_name;
		vector<string> tip_name;
		//vector <Node*> Net_nodes_ptr; /*!< \brief pointers to the nodes of Net \todo use this instead of Net_nodes */
		vector <Node> Net_nodes;  /*!< \brief vector of nodes */
		bool is_net; /*!< \brief true if Net is a network; false if it's a tree */
		bool is_ultrametric; /*!< \brief true if the distances between tips and root are equal; false, otherwise */
		void clear(); 
		void print_all_node();
		void print_all_node_dout();

		Net (){
			string net_str;
			vector <string> tax_name;
			vector <string> tip_name;
			vector <Node> Net_nodes;
			//vector <Node*> Net_nodes_ptr;
		}
		
		//~Net(){
			//tax_name.clear();
			//tip_name.clear();
			//Net_nodes.clear();
		//}
		
		Net(string Net_str);
};

string construct_adding_new_Net_str(Net old_Net);
#endif
