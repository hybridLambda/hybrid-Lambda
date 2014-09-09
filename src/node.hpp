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

#include "utility.hpp"
#include <iomanip>      // std::setw
#include <vector>

#ifndef NODE
#define NODE
/*! \brief Node of a tree or network, it also represent the branch between this node and its parent node
 */

enum NAMETYPE { TAXA, TIP };

class Node {
    friend class Net;
    friend class sim_one_gt;
    friend class HybridLambda;
    friend class Figure;
	public:
        double brchlen1() const { return this->brchlen1_;}
        void set_brchlen1 ( double bl ){ this->brchlen1_ = bl; }
        
        double brchlen2() const { return this->brchlen2_;}
        void set_brchlen2 ( double bl ){ this->brchlen2_ = bl; }

        string label; /*!< \brief String label of a node, each node has unique label */
        size_t node_index; /*!< \brief node index in the array, \todo use this more often!!!*/
        string node_content; /*!< \brief node content, the subtree string at this node */
        bool hybrid; /*!< \brief Hybrid node only, indicator of a hybrid node */
        
    private:    
    //vector<int> descndnt;
	vector<Node*> descndnt_interior_node; /*!< \brief list of pointers to its descndent interior nodes */
	vector<Node*> child; /*!< \brief list of pointers to its child nodes */	
	Node* parent1; /*!< \brief pointer to its parent node. */
	string clade; /*!< \brief clade at this node, \todo this should be modified to a vector <string> */
    
	size_t rank() const { return this->rank_; }
	//int num_child; /*!< \brief number of child \todo this can be replaced by child.size */
	int num_descndnt; /*!< \brief number of the tip nodes, that are descendant from this node */
	int num_descndnt_interior; /*!< \brief number of the interior nodes, that are descendant from this node \todo to be replaced by descndnt_interior_node.size()? */
	vector <double> path_time; 
	double height; /*!< \brief distance to the bottom of the tree */
		
	bool descndnt_of_hybrid; /*!< \brief Indicator of descendant of hybrid nodes. It's true, if it is a descendant of hybrid nodes; false, otherwise. */
	bool tip_bool; /*!< \brief Indicator of tip nodes. It's true, if it is a tip node, otherwise it is false. */
		
	Node* parent2; /*!< \brief Hybrid node only, pointer to its second parent node. */
	//double prob_to_hybrid_left; /*!< \brief Hybrid node only, the probability that a lineage goes to the left */

	string name; /*!< \brief Name of a node, this is not unique for nodes. e.g. if its label is A_1, name is A */
	
	vector <size_t> Net_node_contains_gt_node1; /*!< Used while simulation, check if a Network node contains a gene tree node */
	vector <size_t> Net_node_contains_gt_node2; /*!< Used while simulation, check if a Network node contains a gene tree node */
	
	Node(); /*!< \brief Initialize Node class*/
    void add_child( Node *child_node /*! pointer to the child node*/);
    void CalculateRank();
    
        size_t e_num() const {return this->e_num_;}
        void set_enum( size_t num ) { this->e_num_ = num; }
    
        size_t e_num2() const {return this->e_num2_;}
        void set_enum2( size_t num ) { this->e_num2_ = num; }
    
        bool visited() const { return this->visited_; }
        void set_visited ( bool TorF ){ this->visited_ = TorF; }
    
        void print( bool is_Net );
        void print_dout( bool is_Net );
        void find_tip();
        void find_hybrid_descndnt();
        
        size_t rank_;     /*!< \brief rank of the node, tip node has rank one, the root has the highest rank */        
        bool visited_;
        size_t e_num_;    /*!< \brief numbering the branch */
        size_t e_num2_;   /*!< \brief Hybrid node only, numbering the branch between the node and its second parent */
        double brchlen1_; /*!< \brief Branch length */
        double brchlen2_;/*!< \brief Hybrid node only, Branch length to the second parent*/
        
        bool find_descndnt ( string name, NAMETYPE type );
};


#endif

