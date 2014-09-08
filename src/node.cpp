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

#include"node.hpp"

Node::Node(){
	label=" ";
	node_content=" ";
	num_child=0;
	num_descndnt=0;
	num_descndnt_interior=0;
	parent1=NULL;
	parent2=NULL;
	this->brchlen1_ = 0.0;
	this->brchlen2_ = 0.0;
	this->rank_     = 0;
	this->e_num_    = 0;
	this->e_num2_   = 0;
	//visit=0;
	hybrid=false;
	descndnt_of_hybrid=false;
	tip_bool=false;
	//clade=" ";
	height=0.0;
	//prob_to_hybrid_left=1.0;
	this->visited_ = false;
}


void Node::print( bool is_Net ){
    cout << setw(12) << label;
	if ( is_Net ) cout << setw(6) << hybrid;
    if ( is_Net ) cout << setw(8) << descndnt_of_hybrid;
	cout << setw(5) << tip_bool;
    cout << setw (11) << (this->parent1) ? (parent1->label) : "" ;
	cout << setw (8) << height;
	cout << setw (8) << this->brchlen1();
    if (is_Net){
        cout<<setw (11) << (parent2) ? (parent2->label) : "" ;
        cout<<setw (8) << this->brchlen2();
    }
	cout << setw (7) << num_child;
	cout << setw (8) << num_descndnt;
	cout << setw(4) << num_descndnt_interior;
	cout << setw(6) << this->rank() << "   ";
	//for (size_t i=0;i<descndnt.size();i++){
		//cout<<setw (1)<<descndnt[i];
	//}
	cout << setw(2)<<this->e_num();
	if ( is_Net ) cout << setw(3) << this->e_num2();
	cout << "    " << this->clade;
	//cout<<endl;
}

/*! \brief Add child node to parent node */
void add_node(
Node *parent_node /*! pointer to the parent node*/, 
Node *child_node /*! pointer to the child node*/){
    parent_node->child.push_back(child_node);
	if (child_node->parent1){
		child_node->parent2=parent_node;
		child_node->hybrid=true;
	}
	else {
		child_node->parent1=parent_node;
		//cout<<child_node->label<<" parent1 is "<<child_node->parent1->label<<endl;
	}
	//Node *kidspintr=parent_node->child[parent_node->num_child];
	parent_node->num_child++;
} 


/*! \brief Rank network node from the bottom.
 * 
 * Child node has lower rank than the parent node. Tip nodes have rank one, the root node has the highest rank
 */
void Node::CalculateRank(){
	if ( this->tip_bool ) this->rank_ = 1;
	else {
		size_t child_max_rank = 0;
		for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++ ){
            this->child[ith_child]->CalculateRank();			
			child_max_rank = max( child_max_rank, this->child[ith_child]->rank() );
		}
		this->rank_ = child_max_rank + 1;
	}
}


bool find_descndnt(Node* current, string taxname){
	bool descndnt_found=false;
	if (current->tip_bool){
		if (current->name==taxname){
			descndnt_found=true;
		}
	}
	else {
		for (int i=0;i<current->num_child;i++){
			if (find_descndnt(current->child[i],taxname)){
				descndnt_found=true;
				break;
			}
		}
	}	
	return descndnt_found;
}

bool find_descndnt2(Node* current, string taxname){
	bool descndnt_found=false;
	if (current->tip_bool){
		if (current->label==taxname){
			descndnt_found=true;
		}
		else{
			descndnt_found=false;
		}
	}
	else {
		for (int i=0;i<current->num_child;i++){
			if (find_descndnt2(current->child[i],taxname)){
				descndnt_found=true;
				break;
			}
			else descndnt_found=false;
		}
	}	
	return descndnt_found;
}


/*! \brief label a node if its a tip node */
void Node::find_tip(){
	if ( this->child.size() == 0) this->tip_bool = true;
	else {
		for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++ ){
			(*this->child[ith_child]).find_tip();
            //this->child[ith_child]->find_tip();
		}
	}
}


/*! \brief rewrite node content of nodes */
void rewrite_node_content(vector <Node*> Net_ptr /*! vector of pointers pointing to nodes */){
	int highest_i=0;
	for (size_t i =0; i<Net_ptr.size();i++){
		if (Net_ptr[i]->num_descndnt > Net_ptr[highest_i]->num_descndnt ){highest_i=i;}
	}
	
	//for (size_t node_i=0;node_i<Net_ptr.size();node_i++){
		//Net_ptr[node_i]->print_net_Node();
		//cout<<endl;
	//}
	Net_ptr[highest_i]->CalculateRank();

	for (int rank_i=1;rank_i<=Net_ptr.back()->rank();rank_i++){
		for (size_t node_i=0;node_i<Net_ptr.size();node_i++){
			if (Net_ptr[node_i]->rank()==1){
				Net_ptr[node_i]->node_content=Net_ptr[node_i]->label;
			}
			else{
			if (Net_ptr[node_i]->rank()==rank_i){
				string new_node_content="(";
				for (int child_i=0;child_i<Net_ptr[node_i]->num_child;child_i++){
					ostringstream brchlen_str;
					ostringstream brchlen_str2;
					brchlen_str << Net_ptr[node_i]->child[child_i]->brchlen1();
					if (Net_ptr[node_i]->child[child_i]->node_content==Net_ptr[node_i]->child[child_i]->label){
						new_node_content=new_node_content+Net_ptr[node_i]->child[child_i]->label+":"+brchlen_str.str();}
					else{
						bool new_hybrid_node=false;
						for (size_t node_ii=0;node_ii<node_i;node_ii++){
							for (int node_ii_child_i=0;node_ii_child_i<Net_ptr[node_ii]->num_child;node_ii_child_i++){
								if (Net_ptr[node_ii]->child[node_ii_child_i]->node_content==Net_ptr[node_i]->child[child_i]->node_content){
									new_hybrid_node=true;
									brchlen_str2 << Net_ptr[node_i]->child[child_i]->brchlen2();
								break;}
							}
							//if (new_hybrid_node==1){break;}
						}
						if (new_hybrid_node){
							new_node_content=new_node_content+Net_ptr[node_i]->child[child_i]->label+":"+brchlen_str2.str();
						}
						else{
							new_node_content=new_node_content+Net_ptr[node_i]->child[child_i]->node_content+Net_ptr[node_i]->child[child_i]->label+":"+brchlen_str.str();
						}
					}
					//new_node_content=new_node_content+sp_nodes_ptr_rm1[node_i]->child[child_i]->node_content+sp_nodes_ptr_rm1[node_i]->child[child_i]->label+":"+brchlen_str.str();
					if (child_i<Net_ptr[node_i]->num_child-1){
						new_node_content=new_node_content+",";
					}
				}
				new_node_content=new_node_content+")";
				Net_ptr[node_i]->node_content=new_node_content;
				}
			}
		}	
	}
	//for (size_t i =0; i<Net_ptr.size();i++){
		//cout<<Net_ptr[i]->label<<" "<<Net_ptr[i]->node_content<<endl;//<<"!!!!";
	//}
}


///////////////////////////////////////////////////////////////////////////////////////////// consider removed

/*! \brief Label a node if its a descendant of a hybrid node */
void Node::find_hybrid_descndnt(){
	if ( this->tip_bool ) return;
    for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++){
        if ( this->hybrid || this->descndnt_of_hybrid ) this->child[ith_child]->descndnt_of_hybrid = true;
        this->child[ith_child]->find_hybrid_descndnt();
    }
}
