//node.cpp
#include"node.hpp"



Node::Node(){
	label=" ";
	node_content=" ";
	num_child=0;
	num_descndnt=0;
	num_descndnt_interior=0;
	parent1=NULL;
	parent2=NULL;
	brchlen1=0.0;
	brchlen2=0.0;
	rank=0;
	e_num=0;
	e_num2=0;
	//visit=0;
	hybrid=false;
	e_num2=0;
	descndnt_of_hybrid=false;
	tip_bool=false;
	//clade=" ";
	absolute_time=0.0;
	//prob_to_hybrid_left=1.0;
	visited = false;
	//double pAC=0.0;
}




void Node::print_net_Node(){
	cout<<setw(12)<<label;
	cout<<setw(6)<<hybrid;
	cout<<setw(8)<<descndnt_of_hybrid;
	cout<<setw(5)<<tip_bool;
	if (parent1){cout<<setw (11)<<(parent1->label);}
	else cout<<"           ";
	cout<<setw (8)<<absolute_time;
	cout<<setw (8)<<brchlen1;
	if (parent2){cout<<setw (9)<<parent2->label;}
	else cout<<"         ";
	cout<<setw (8)<<brchlen2;
	cout<<setw (7)<<num_child;
	cout<<setw (8)<<num_descndnt;
	cout<<setw(4)<<num_descndnt_interior;
//	cout<<setw (7)<<current.e_num;
//	cout<<setw (3)<<current.e_num2;
	cout<<setw (6)<<rank<<"   ";
	//for (unsigned int i=0;i<descndnt.size();i++){
		//cout<<setw (1)<<descndnt[i];
	//}
	cout<<setw(2)<<e_num;
	cout<<setw(3)<<e_num2;
	cout<<"    "<<clade;
	//cout<<endl;
}


void Node::print_tree_Node(){
	cout<<setw (12)<<label;
	cout<<setw(5)<<tip_bool;
	if (parent1){cout<<setw (11)<<(parent1->label);}
	else cout<<"           ";
	cout<<setw (11)<<absolute_time;
	cout<<setw (11)<<brchlen1;
	cout<<setw (7)<<num_child;
	cout<<setw (8)<<num_descndnt;
	cout<<setw(4)<<num_descndnt_interior;
	cout<<setw (6)<<rank<<"   ";
	//for (unsigned int i=0;i<descndnt.size();i++){
		//cout<<setw (1)<<descndnt[i];
	//}	
	cout<<setw(3)<<e_num;
	cout<<"    "<<clade;
	//cout<<setw(3)<<num_descndnt_interior;
}

	//~Node(){
void Node::clear(){
	//clade=" ";
	clade.clear();
	//~clade();
	label=" ";
	node_content=" ";
	num_child=0;
	num_descndnt=0;
	parent1=NULL;
	parent2=NULL;
	brchlen1=0.0;
	brchlen2=0.0;
	rank=0;
	e_num=0;
	//visit=0;
	hybrid=false;
	e_num2=0;
	descndnt_of_hybrid=false;
	tip_bool=false;
	label.clear();
	node_content.clear();
	child.clear();
	//descndnt.clear();
}


void print_all_child(Node *parent /*! pointer to the parent node*/){
    cout << parent->label << " has " << parent->num_child << " kids" << endl;
    for (int i_num_child=0;i_num_child<=parent->num_child-1;i_num_child++){
        cout<<parent->child[i_num_child]->label<<endl;
    }
}


void print_parent(Node *child /*! pointer to the child node*/){
	cout<<child->label<<" has parents "<<endl;
	cout<<child->parent1->label<<" and "<<child->parent2->label <<endl;
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
int ranking(Node *current){
	int current_rank;
	if (current->tip_bool){
		current->rank=1;}
	else
	{
		int child_max_rank=0;
		for (int i_num_child=0;i_num_child<current->num_child;i_num_child++){
			int child_rank=ranking(current->child[i_num_child]);
			child_max_rank=max(child_rank,child_max_rank);
		}
		current->rank=child_max_rank+1;
	}
		
	//cout<<current->label<<"  "<<current->rank<<endl;
	
	return current_rank=current->rank;
}


bool find_descndnt(Node* current, string taxname){
	bool descndnt_found=false;
	if (current->tip_bool){
		//if (current->label==taxname){
			//cout<<current->name<<endl;
			//cout<<taxname<<endl;
		if (current->name==taxname){
			descndnt_found=true;
		}
		//else{
			//descndnt_found=false;
		//}
	}
	else {//int i;
		for (int i=0;i<current->num_child;i++){
			if (find_descndnt(current->child[i],taxname)){
				descndnt_found=true;
				break;
			}
			//else descndnt_found=false;
		}
	}	
	return descndnt_found;
}

bool find_descndnt2(Node* current, string taxname){
	bool descndnt_found=false;
	if (current->tip_bool){
		//if (current->label==taxname){
		if (current->label==taxname){
			descndnt_found=true;
		}
		else{
			descndnt_found=false;
		}
	}
	else {//int i;
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
void find_tip(Node *current /*! pointer to a node*/){
	if (current->num_child==0){
	//cout<<current->label<<" is a tip "<<endl;
		current->tip_bool=true;
	}
	else {
		for (int i_num_child=0;i_num_child < current->num_child;i_num_child++){
			//if (current->hybrid || current->descndnt_of_hybrid){
				//current->child[i_num_child]->descndnt_of_hybrid=true;
			//}
			find_tip(current->child[i_num_child]);
		}
	}
}


/*! \brief Label a node if its a descendant of a hybrid node */
void find_hybrid_descndnt(Node *current/*! pointer to a node*/){
	//if (current->num_child==0){
		//current->tip_bool=false;}
	//else {
	if (!current->tip_bool){
		for (int i_num_child=0;i_num_child<current->num_child;i_num_child++){
			if (current->hybrid || current->descndnt_of_hybrid){
				current->child[i_num_child]->descndnt_of_hybrid=true;
			}
			find_hybrid_descndnt(current->child[i_num_child]);
		}
	}
}




/*! \brief rewrite node content of nodes */
void rewrite_node_content(vector <Node*> Net_ptr /*! vector of pointers pointing to nodes */){
	int highest_i=0;
	for (unsigned int i =0; i<Net_ptr.size();i++){
		if (Net_ptr[i]->num_descndnt > Net_ptr[highest_i]->num_descndnt ){highest_i=i;}
	}
	
	//for (unsigned int node_i=0;node_i<Net_ptr.size();node_i++){
		//Net_ptr[node_i]->print_net_Node();
		//cout<<endl;
	//}
	
	ranking(Net_ptr[highest_i]);
	//cout<<Net_ptr[highest_i]->node_content<<endl;
	//for (unsigned int node_i=0;node_i<Net_ptr.size();node_i++){
		//Net_ptr[node_i]->print_net_Node();
		//cout<<endl;
	//}
	for (int rank_i=1;rank_i<=Net_ptr.back()->rank;rank_i++){
		for (unsigned int node_i=0;node_i<Net_ptr.size();node_i++){
			if (Net_ptr[node_i]->rank==1){
				Net_ptr[node_i]->node_content=Net_ptr[node_i]->label;
			}
			else{
			if (Net_ptr[node_i]->rank==rank_i){
				string new_node_content="(";
				for (int child_i=0;child_i<Net_ptr[node_i]->num_child;child_i++){
					ostringstream brchlen_str;
					ostringstream brchlen_str2;
					brchlen_str<<Net_ptr[node_i]->child[child_i]->brchlen1;
					if (Net_ptr[node_i]->child[child_i]->node_content==Net_ptr[node_i]->child[child_i]->label){
						new_node_content=new_node_content+Net_ptr[node_i]->child[child_i]->label+":"+brchlen_str.str();}
					else{
						bool new_hybrid_node=false;
						for (unsigned int node_ii=0;node_ii<node_i;node_ii++){
							for (int node_ii_child_i=0;node_ii_child_i<Net_ptr[node_ii]->num_child;node_ii_child_i++){
								if (Net_ptr[node_ii]->child[node_ii_child_i]->node_content==Net_ptr[node_i]->child[child_i]->node_content){
									new_hybrid_node=true;
									brchlen_str2<<Net_ptr[node_i]->child[child_i]->brchlen2;
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
	//for (unsigned int i =0; i<Net_ptr.size();i++){
		//cout<<Net_ptr[i]->label<<" "<<Net_ptr[i]->node_content<<endl;//<<"!!!!";
	//}
}

