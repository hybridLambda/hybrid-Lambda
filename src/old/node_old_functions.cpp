
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


void Node::print_net_Node(){
	cout<<setw(12)<<label;
	cout<<setw(6)<<hybrid;
	cout<<setw(8)<<descndnt_of_hybrid;
	cout<<setw(5)<<tip_bool;
	//if (parent1){cout<<setw (11)<<(parent1->label);}
	//else cout<<"           ";
    cout<<setw (11)<< (this->parent1) ? (parent1->label) : "" ;
	cout<<setw (8)<<height;
	cout<<setw (8)<<brchlen1;
	//if (parent2){cout<<setw (9)<<parent2->label;}
	//else cout<<"         ";
    cout<<setw (11)<< (parent2) ? (parent2->label) : "" ;

	cout<<setw (8)<<brchlen2;
	cout<<setw (7)<<num_child;
	cout<<setw (8)<<num_descndnt;
	cout<<setw(4)<<num_descndnt_interior;
//	cout<<setw (7)<<current.e_num;
//	cout<<setw (3)<<current.e_num2;
	cout<<setw (6)<<rank<<"   ";
	//for (size_t i=0;i<descndnt.size();i++){
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
    cout<<setw (11)<< (parent1) ? (parent1->label) : "" ;
	cout<<setw (11)<<height;
	cout<<setw (11)<<brchlen1;
	cout<<setw (7)<<num_child;
	cout<<setw (8)<<num_descndnt;
	cout<<setw(4)<<num_descndnt_interior;
	cout<<setw (6)<<rank<<"   ";
	//for (size_t i=0;i<descndnt.size();i++){
		//cout<<setw (1)<<descndnt[i];
	//}	
	cout<<setw(3)<<e_num;
	cout<<"    "<<clade;
	//cout<<setw(3)<<num_descndnt_interior;
}


int ranking(Node *current){
	int current_rank;
	if (current->tip_bool){
		current->rank=1;}
	else
	{
		int child_max_rank=0;
		for (int ith_child=0;ith_child<current->num_child;ith_child++){
			int child_rank=ranking(current->child[ith_child]);
			child_max_rank=max(child_rank,child_max_rank);
		}
		current->rank=child_max_rank+1;
	}
		
	//cout<<current->label<<"  "<<current->rank<<endl;
	
	return current_rank=current->rank;
}
