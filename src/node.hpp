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
//node.hpp
#include"utility.hpp"
#ifndef NODE
#define NODE
/*! \brief Node of a tree or network, it also represent the branch between this node and its parent node
 */

class Node {
	public:
	//vector<int> descndnt;
	vector<Node*> descndnt_interior_node; /*!< \brief list of pointers to its descndent interior nodes */
	vector<Node*> child; /*!< \brief list of pointers to its child nodes */	
	Node* parent1; /*!< \brief pointer to its parent node. */
	string clade; /*!< \brief clade at this node, \todo this should be modified to a vector <string> */
	string label; /*!< \brief String label of a node, each node has unique label */
	string node_content; /*!< \brief node content, the subtree string at this node */
//	class Net * SubNetwork; /*!< \brief \todo pointer to the subtree of this node */
	unsigned int e_num; /*!< \brief numbering the branch */
	int rank; /*!< \brief rank of the node, tip node has rank one, the root has the highest rank */
	int num_child; /*!< \brief number of child \todo this can be replaced by child.size */
	int num_descndnt; /*!< \brief number of the tip nodes, that are descendant from this node */
	int num_descndnt_interior; /*!< \brief number of the interior nodes, that are descendant from this node \todo to be replaced by descndnt_interior_node.size()? */
	vector <double> path_time; 
	double absolute_time; /*!< \brief distance to the bottom of the tree */
	double brchlen1; /*!< \brief Branch length */
	bool visited;
	bool descndnt_of_hybrid; /*!< \brief Indicator of descendant of hybrid nodes. It's true, if it is a descendant of hybrid nodes; false, otherwise. */
	bool tip_bool; /*!< \brief Indicator of tip nodes. It's true, if it is a tip node, otherwise it is false. */
	
	unsigned int node_index; /*!< \brief node index in the array, \todo use this more often!!!*/
	
	
	/* These members apply to only hybrid nodes */
	bool hybrid; /*!< \brief Hybrid node only, indicator of a hybrid node */
	Node* parent2; /*!< \brief Hybrid node only, pointer to its second parent node. */
	unsigned int e_num2; /*!< \brief Hybrid node only, numbering the branch between the node and its second parent */
	//double prob_to_hybrid_left; /*!< \brief Hybrid node only, the probability that a lineage goes to the left */
	double brchlen2;/*!< \brief Hybrid node only, Branch length to the second parent*/

	string name; /*!< \brief Name of a node, this is not unique for nodes. e.g. if its label is A_1, name is A */
	
	vector <unsigned int> Net_node_contains_gt_node1; /*!< Used while simulation, check if a Network node contains a gene tree node */
	vector <unsigned int> Net_node_contains_gt_node2; /*!< Used while simulation, check if a Network node contains a gene tree node */
	
	Node(); /*!< \brief Initialize Node class*/
	void print_net_Node();
	void print_tree_Node();
		void print_tree_Node_dout();
void print_net_Node_dout();
	void clear();
	
		
};

void add_node(Node *parent_node, Node *child_node);
void find_tip(Node *current);
void find_hybrid_descndnt(Node *current);
bool find_descndnt(Node* current, string taxname);
bool find_descndnt2(Node* current, string taxname);
void rewrite_node_content(vector <Node*> Net_ptr);
int ranking(Node *current);

#endif

