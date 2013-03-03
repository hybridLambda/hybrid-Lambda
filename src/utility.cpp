/* 
 * 
 * hybrid_sim is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_sim 
 * 
 * hybrid_sim is free software: you can redistribute it and/or modify
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

/*! \file utility.cpp
 *  \brief Core function of converting a Newick (extended Newick) format string into a species tree (network), and simple string manipulation for tree strings
 */



#include"utility.hpp"
extern bool debug_bool;
//extern bool log_file_bool;


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

/*! \brief Identify if its the start of the taxon name in a newick string, should be replaced by using (isalpha() || isdigit())  */
bool start_of_tax_name(string in_str,size_t i){
	bool start_bool=false;
	if ( (in_str[i]!='(' && in_str[i-1]=='(') || (in_str[i-1]==',' && in_str[i]!='(') || ( (in_str[i-1]==')') && ( in_str[i]!=')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) ) {
		start_bool=true;	
	}
	
	return 	start_bool;
}

/*! \brief Checking Parenthesis of a (extended) Newick string */
void checking_Parenthesis(string in_str){
	int num_b=0;
	for (size_t i=0;i<in_str.size();i++){
		if (in_str[i]=='('){
			num_b++;
		}
		if (in_str[i]==')'){
			num_b--;
		}
	}
	if (num_b!=0){
		cout<<"Error:"<<endl;
		cout<<in_str<<endl;
		cout<<"Parenthesis not balanced!"<<endl;
		exit (1);
	}
}

size_t Parenthesis_balance_index_backwards(string in_str,size_t i){
	size_t j=i;
	int num_b=0;
	for (;j>0;j--){
		if (in_str[j]=='('){
			num_b--;
		}
		if (in_str[j]==')'){
			num_b++;
		}
		if (num_b==0){
			break;
		}
	}
	return j;
}

size_t Parenthesis_balance_index_forwards(string in_str,size_t i){
	size_t j=i;
	int num_b=0;
	for (;j<in_str.size();j++){
		if (in_str[j]=='('){
			num_b++;
		}
		if (in_str[j]==')'){
			num_b--;
		}
		if (num_b==0){
			break;
		}
	}
	return j;
}



/*! \brief Construct Net object from a (extended) Newick string */
Net::Net(string old_string /*! input (extended) newick form string */){
	if (old_string.size()>0){

		checking_Parenthesis(old_string);
		net_str=checking_labeled(old_string);
	
		if (debug_bool){
			cout<<"Net::Net flag1"<<endl;
		}

		int net_str_len=net_str.length();
		vector<string> labels;
		vector<string> node_contents;
		vector<string> brchlens;
		size_t found_bl=net_str.find(':');
		for (size_t i_str_len=1;i_str_len<net_str.size();){
			if (net_str[i_str_len]=='e' && (net_str[i_str_len+1]=='-' || net_str[i_str_len+1]=='+')){
				i_str_len++;
			}
			else{
				//if (isalpha(net_str[i_str_len])){
				if ( start_of_tax_name(net_str,i_str_len)){
				//if (isalpha(net_str[i_str_len]) || isdigit(net_str[i_str_len])){	
					size_t str_start_index=i_str_len;
					string label=extract_label(net_str,i_str_len);
					//cout<<label<<endl;
					labels.push_back(label);
					//int str_end_index=label.size()+i_str_len-1;
					string node_content;
					if (net_str[str_start_index-1]==')'){
						size_t rev_dummy_i=Parenthesis_balance_index_backwards(net_str,str_start_index-1);
						size_t substr_len=str_start_index-rev_dummy_i;
						node_content=net_str.substr(rev_dummy_i,substr_len);			
					}
					else {
						node_content=label;
					}
					i_str_len=label.size()+i_str_len;
					//cout<<node_content<<endl;
					node_contents.push_back(node_content);
					//size_t found_bl=net_str.find(':');
					string brchlen;
					if (found_bl!=string::npos){					
						//for (size_t i_num_str=i_str_len+1;i_num_str<net_str_len;i_num_str++){	
							//brchlen.push_back(net_str[i_num_str]);			
							//if (net_str[i_num_str+1]==',' || net_str[i_num_str+1]==')'){
								//break;}
								
						//}
						size_t found=min(min(net_str.find(",",i_str_len+1),net_str.find(")",i_str_len+1)),net_str.size());
						brchlen=net_str.substr(i_str_len+1,found-i_str_len-1);
					
					}
					found_bl=net_str.find(":",found_bl+1);
					brchlens.push_back(brchlen);
				}
				else {
					i_str_len++;
				}
			}
		}
			
		if (debug_bool){
			cout<<"Net::Net flag2"<<endl;
		}
		//vector <Node> Net_nodes;
		int label_counter=brchlens.size();
		for (int new_i_label=0;new_i_label<label_counter;new_i_label++){
			Node empty_node;
			//empty_node.pAC=0.0;
			Net_nodes.push_back(empty_node);
			Net_nodes[new_i_label].label=labels[new_i_label];
			Net_nodes[new_i_label].node_content=node_contents[new_i_label];
			//cout<<Net_nodes[new_i_label].label<<" "<<Net_nodes[new_i_label].node_content<<endl;
			string s(brchlens[new_i_label]);
			istringstream istr(s);
			istr>>Net_nodes[new_i_label].brchlen1;
		}
		//int repeated_num_node=0;
		for (unsigned int i=1;i<Net_nodes.size()-1;i++){
			unsigned int j;
			for ( j=i+1;j<Net_nodes.size()-1;j++){
				if (Net_nodes[j].label==Net_nodes[i].label){
					//repeated_num_node++;
					if (Net_nodes[j].node_content[0]=='('){
						Net_nodes[i].node_content=Net_nodes[j].node_content;
	//					Net_nodes[i].brchlen2=Net_nodes[i].brchlen1;
	//					double* brch_ptr_i=&Net_nodes[i].brchlen1;
	//					double* brch_ptr_j=&Net_nodes[j].brchlen1;
	//					brch_ptr_i=brch_ptr_j;
					}
	//				else{
					Net_nodes[i].brchlen2=Net_nodes[j].brchlen1;
	//				}
					break;
				}
			}
			if (Net_nodes[j].label==Net_nodes[i].label){
				Net_nodes.erase(Net_nodes.begin()+j);
			//	Net_nodes_ptr.erase(Net_nodes_ptr.begin()+j);
			}
		}
		//bool multi_label_bool=false;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			if(Net_nodes[i].label==Net_nodes[i].node_content){
				if (Net_nodes[i].label.find("_")>0){
					//multi_label_bool=true;
					Net_nodes[i].name=Net_nodes[i].label.substr(0,Net_nodes[i].label.find("_"));
					//cout<<Net_nodes[i].name<<endl;
					bool new_tax_bool=true;
					for (unsigned int tax_i=0;tax_i<tax_name.size();tax_i++){
						if (tax_name[tax_i]==Net_nodes[i].name){
							new_tax_bool=false;
							break;
						}
					}
					if (new_tax_bool){
						tax_name.push_back(Net_nodes[i].name);
					}
					//cout<<tax_name.back()<<endl;
				}
				else{
					tax_name.push_back(Net_nodes[i].label);
				}
				tip_name.push_back(Net_nodes[i].label);
			}
		}
		sort(tax_name.begin(), tax_name.end());
		sort(tip_name.begin(), tip_name.end());
		//for (unsigned int i=0;i<Net_nodes.size();i++){
				//valarray <int> intial_descndnt(0,tax_name.size());
				//descndnt.push_back(intial_descndnt);
		//}
		if (debug_bool){
			cout<<"Net::Net flag3"<<endl;
		}
		vector <Node*> Net_nodes_ptr;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			Net_nodes[i].node_index=i;
			Node* new_node_ptr=NULL;
			Net_nodes_ptr.push_back(new_node_ptr);
			Net_nodes_ptr[i]=&Net_nodes[i];
		}
		for (unsigned int i=0;i<Net_nodes.size();i++){
			if (Net_nodes[i].node_content[0]=='('){
				char child_node1[Net_nodes[i].node_content.length()];
				unsigned int i_content_len;
				unsigned int j_content_len;
				for (i_content_len=1;i_content_len<Net_nodes[i].node_content.length();){
					//if (Net_nodes[i].node_content[i_content_len]=='(' ||  isalpha(Net_nodes[i].node_content[i_content_len])){	
					if (Net_nodes[i].node_content[i_content_len]=='(' ||  start_of_tax_name(Net_nodes[i].node_content,i_content_len) ){	
					//if (Net_nodes[i].node_content[i_content_len]=='(' ||  isalpha(Net_nodes[i].node_content[i_content_len]) || isdigit(Net_nodes[i].node_content[i_content_len]) ){	
						if (Net_nodes[i].node_content[i_content_len]=='('){
							j_content_len=Parenthesis_balance_index_forwards(Net_nodes[i].node_content,i_content_len)+1;
						}
						else{
							j_content_len=i_content_len;
						}
						int child1_node_content_i=0;
						for (;j_content_len<Net_nodes[i].node_content.length(); j_content_len++){
							child_node1[child1_node_content_i]=Net_nodes[i].node_content[j_content_len];
							char stop=Net_nodes[i].node_content[j_content_len+1];
							if (stop==',' || stop==')' || stop==':'){
								child_node1[child1_node_content_i+1]='\0';
								break;}
							child1_node_content_i++;
							}
							string child_node1_str=child_node1;		
							i_content_len=j_content_len+2;
							for (unsigned int j=0;j<Net_nodes.size();j++){
								if (child_node1_str==Net_nodes[j].label){
									add_node(Net_nodes_ptr[i],Net_nodes_ptr[j]);
								}
							}
					}
					else{i_content_len++;}
				}	
			}
		}
		//cout<<descndnt.size()<<endl;
		
		if (debug_bool){
			cout<<"Net::Net flag4"<<endl;
		}
		find_tip(Net_nodes_ptr.back());
		//cout<<"number of child " <<Net_nodes_ptr.back()->child.size()<<endl;
		//cout<<"node content "<<Net_nodes_ptr.back()->node_content<<endl;
		find_hybrid_descndnt(Net_nodes_ptr.back());
		max_rank=ranking(Net_nodes_ptr.back());
	
		for (unsigned int i=0;i<Net_nodes.size();i++){
			valarray <int> descndnt_dummy(0,tax_name.size());
			descndnt.push_back(descndnt_dummy);
			valarray <int> descndnt2_dummy(0,tip_name.size());
			descndnt2.push_back(descndnt2_dummy);
			//cout<<Net_nodes_ptr[i]->name<<"  "<<Net_nodes_ptr[i]->label<<"  "<<Net_nodes_ptr[i]->tip_bool << " ";
			for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				//cout<<Net_nodes_ptr[i]->label<<"   "<<tax_name[tax_name_i]<<"  ";
				if (find_descndnt(Net_nodes_ptr[i],tax_name[tax_name_i])){
					descndnt[i][tax_name_i]=1;
				}
				//cout<<descndnt[i][tax_name_i];
			}
			//cout<<endl;
			for (unsigned int tip_name_i=0;tip_name_i<tip_name.size();tip_name_i++){
				if (find_descndnt2(Net_nodes_ptr[i],tip_name[tip_name_i])){
					descndnt2[i][tip_name_i]=1;
				}
			}
			Net_nodes[i].num_descndnt=descndnt[i].sum();
		}
		if (debug_bool){
			cout<<"Net::Net flag5"<<endl;
		}
		
		//num_descndnt_interior(Net_nodes_ptr.back()); //replace by following
		for (unsigned int i=0;i<Net_nodes_ptr.size();i++){
			for (unsigned int j=0;j<Net_nodes_ptr.size();j++){
				if (i!=j){
					valarray <int> descndnt_diff=(descndnt[i]-descndnt[j]);
					//cout<<"hea"<<endl;
					if (descndnt_diff.min() >= 0 && Net_nodes[i].rank > Net_nodes[j].rank && Net_nodes[j].rank>=2){
						Net_nodes_ptr[i]->num_descndnt_interior=Net_nodes_ptr[i]->num_descndnt_interior+1;
						Net_nodes_ptr[i]->descndnt_interior_node.push_back(Net_nodes_ptr[j]);
					}
					
				}
			
			}
			//cout<<"checking #interior_des "<<Net_nodes_ptr[i]->num_descndnt_interior<<" "<<Net_nodes_ptr[i]->descndnt_interior_node.size()<<endl;
		}
		if (debug_bool){
			cout<<"Net::Net flag6"<<endl;
		}	
		int e_num_old=0;
		enumerate_internal_branch(Net_nodes_ptr.back(),e_num_old);
		if (debug_bool){
			cout<<"Net::Net flag6.5"<<endl;
		}
		//cout<<descndnt.size()<<endl;
		//cout<<net_str<<endl;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			//Net_nodes[i].print_tree_Node();
			//cout<<endl;
			//for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				//cout<<descndnt[i][tax_name_i];
			//}
			//cout<<endl;
			//cout<<i<<" "<<Net_nodes[i].name<<" " <<Net_nodes[i].label<<Net_nodes[i].clade<<endl;
			//cout<<tax_name.size()<<endl;
			for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				//cout<<tax_name[tax_name_i]<<endl;
				if (descndnt[i][tax_name_i] == 1){
					//cout<<i<<" "<<Net_nodes[i].label<<" "<<Net_nodes[i].clade<<endl;
					if (Net_nodes[i].clade.size()==0){
						Net_nodes[i].clade=tax_name[tax_name_i];
					}
					else{
						Net_nodes[i].clade=Net_nodes[i].clade+tax_name[tax_name_i];
					}
					Net_nodes[i].clade.push_back('&');
				}
			}
			//cout<<i<<" "<<Net_nodes[i].clade<<endl;
			
			Net_nodes[i].clade.erase(Net_nodes[i].clade.size()-1,1);
		}
		
		if (debug_bool){
			cout<<"Net::Net flag7"<<endl;
		}
		//check for coaleased tips(& sign in the tips)
		bool rewrite_descndnt=false;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			if (Net_nodes[i].tip_bool ){
				for (unsigned int i_str=0;i_str<Net_nodes[i].clade.size();i_str++){
					if (Net_nodes[i].clade[i_str]=='&'){
						rewrite_descndnt=true;
						break;
					}
				}
			}
			if (rewrite_descndnt){
				break;
			}
		}
		
		if (rewrite_descndnt){
			tax_name.clear();
			int tax_name_start=0;
			int tax_name_length=0;
			for (unsigned int new_i_str=0;new_i_str<Net_nodes.back().clade.size();new_i_str++){
				tax_name_length++;
				if (Net_nodes.back().clade[new_i_str]=='&'){
					tax_name_length--;
					tax_name.push_back(Net_nodes.back().clade.substr(tax_name_start,tax_name_length));
					tax_name_start=new_i_str+1;
					tax_name_length=0;
				}				
				if (new_i_str==Net_nodes.back().clade.size()-1){
					tax_name.push_back(Net_nodes.back().clade.substr(tax_name_start,tax_name_length));
				}
			}
			sort(tax_name.begin(), tax_name.end());
			//cout<<descndnt.size()<<endl;
			descndnt.clear();
		//	cout<<descndnt.size()<<endl;
			//~descndnt();
			for (unsigned int i=0;i<Net_nodes.size();i++){
				vector <string> contained_tips;
				valarray <int> re_initial_descndnt(0,tax_name.size());
				int tax_name_start=0;
				int tax_name_length=0;
				for (unsigned int new_i_str=0;new_i_str<Net_nodes[i].clade.size();new_i_str++){
					tax_name_length++;
					if (Net_nodes.back().clade[new_i_str]=='&'){
						tax_name_length--;
						contained_tips.push_back(Net_nodes[i].clade.substr(tax_name_start,tax_name_length));
						tax_name_start=new_i_str+1;
						tax_name_length=0;
					}				
					if (new_i_str==Net_nodes[i].clade.size()-1){
						contained_tips.push_back(Net_nodes[i].clade.substr(tax_name_start,tax_name_length));
					}
				}
				for (unsigned int tax_i=0;tax_i<tax_name.size();tax_i++){
					for (unsigned int contained_tax_i=0;contained_tax_i<contained_tips.size();contained_tax_i++){
						if (tax_name[tax_i]==contained_tips[contained_tax_i]){
							//descndnt[i][tax_i]=1;
							re_initial_descndnt[tax_i]=1;
						}
					}
				}	
				descndnt.push_back(re_initial_descndnt);
			}
			
			for (unsigned int i=0;i<Net_nodes.size();i++){
				Net_nodes[i].clade=" ";
				for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
					if (descndnt[i][tax_name_i] == 1){
						if (Net_nodes[i].clade == " "){
							Net_nodes[i].clade=tax_name[tax_name_i];
						}
						else{
							Net_nodes[i].clade=Net_nodes[i].clade+tax_name[tax_name_i];
						}
						Net_nodes[i].clade.push_back('&');
					}
				}
				Net_nodes[i].clade.erase(Net_nodes[i].clade.size()-1,1);				
			}
		}
		
		is_net=is_net_func();
		//is_ultrametric=true;
		is_ultrametric=is_ultrametric_func();

	}
	else{
		descndnt.clear();
		tax_name.clear();
		Net_nodes.clear();
	}
		
	if (debug_bool){
		cout<<"Net::Net flag8"<<endl;
	}
}

string Net::checking_labeled(string in_str){
	bool labeled_bool=true;
	string out_str;
	for (size_t i=0;i<in_str.size();i++){
		if (in_str[i]==')' && i==end_of_label_or_bl(in_str, i) ){
			labeled_bool=false;
			break;
		}
	}
	
	if (labeled_bool){
		out_str=in_str;
	}
	else{
		out_str=label_interior_node(in_str);
	}	
	return out_str;
}


bool Net::is_ultrametric_func(){
	bool is_ultrametric_return=true;
	vector <int> remaining_node(Net_nodes.size(),0);
	for (unsigned int node_i=0;node_i<Net_nodes.size();node_i++){
		remaining_node[node_i]=node_i;
	
	}
	//cout<<remaining_node.size()<<endl;
	//for (int rank_i=1;rank_i<=max_rank;rank_i++){
		//for (unsigned int node_i=0;node_i<Net_nodes.size();node_i++){
	int rank_i=1;
	unsigned int remaining_node_i=0;	
	while (remaining_node.size()>0){
		//remaining_node_i=0;
		int node_i=remaining_node[remaining_node_i];
		//cout<<remaining_node[remaining_node_i]<<endl;
		//rank_i=1;
		if (Net_nodes[node_i].rank==rank_i){
			if (rank_i==1){
				Net_nodes[node_i].path_time.push_back(0.0);
			}
			else{
				for (unsigned int child_i=0;child_i<Net_nodes[node_i].child.size();child_i++){
					double current_child_time;
					if (Net_nodes[node_i].child[child_i]->parent1->label==Net_nodes[node_i].label){
						current_child_time=Net_nodes[node_i].child[child_i]->brchlen1;
					}
					else{
						//if (Net_nodes[node_i].child[child_i]->parent2->label==Net_nodes[node_i].label){
							current_child_time=Net_nodes[node_i].child[child_i]->brchlen2;
						//}
						//else{
							//cout<<"warning!!!!! check code again"<<endl;
						//}
							
					}
					for (unsigned int child_i_time_i=0;child_i_time_i<Net_nodes[node_i].child[child_i]->path_time.size();child_i_time_i++){
						Net_nodes[node_i].path_time.push_back(current_child_time+Net_nodes[node_i].child[child_i]->path_time[child_i_time_i]);
					}
				}
			}
			
			remaining_node.erase(remaining_node.begin()+remaining_node_i);
			//cout<<rank_i<<" "<<Net_nodes[node_i].label<<" "<<Net_nodes[node_i].path_time.size()<<endl;
		}
		else{
			remaining_node_i++;
		}
		//
		if (remaining_node_i==remaining_node.size()-1){
			rank_i++;
			remaining_node_i=0;
		}
	}
	//}
	for (unsigned int node_i=0;node_i<Net_nodes.size();node_i++){
		//bool node_i_path_time_eq=true;
		//cout<<Net_nodes[node_i].label<<"  ";
		for (unsigned int path_time_i=0;path_time_i<Net_nodes[node_i].path_time.size();path_time_i++){
			if (pow((Net_nodes[node_i].path_time[path_time_i]-Net_nodes[node_i].path_time[0]),2)>0.000001){
				is_ultrametric_return=false;
				//cout<<Net_nodes[node_i].label<<endl;
				break;
			}
			//cout<<Net_nodes[node_i].path_time[path_time_i]<<" ";
		}
		//cout<<Net_nodes[node_i].label<<endl;
		//if (!is_ultrametric){
			//break;
		//}
		Net_nodes[node_i].absolute_time=Net_nodes[node_i].path_time[0];
	}
	return is_ultrametric_return;
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


//void label_des_indictr(Node *parent, int des_indictr[], string tax_name[], int tax_counter){
    //if (parent->num_child==0){
		//for (int i=0;i<tax_counter;i++){
			//if (parent->label==tax_name[i] && des_indictr[i]==0){
				//des_indictr[i]=des_indictr[i]+1;
			//}
		//}
	//}
	//for (int i_num_child=0;i_num_child<parent->num_child;i_num_child++){
		//label_des_indictr(parent->child[i_num_child], des_indictr, tax_name, tax_counter);
    //}
//}


//int compute_num_descndnt_interior(Node* current){
	//int myreturn=0;
	//if (current->rank <= 2){
		//current->num_descndnt_interior=0;
	//}
	//else{
		////int total_num_descndnt_interior=0;
		//for (int i=0;i<current->num_child;i++){
			//if (current->child[i]->rank==2){
				//current->num_descndnt_interior=current->num_descndnt_interior+1;
			//}
			//else{
				//if (current->child[i]->rank>2){
					//current->num_descndnt_interior=1+current->num_descndnt_interior+num_descndnt_interior(current->child[i]);
				//}
				////current->child[i]->num_descndnt_interior=current->child[i]->num_descndnt_interior+1;
				////child_num_descndnt_interior=num_descndnt_interior(current->child[i_num_child]);
			//}
			
			////total_num_descndnt_interior=total_num_descndnt_interior+child_num_descndnt_interior;
		//}
		////current->num_descndnt_interior=total_num_descndnt_interior	
	//}
	
	
	//return myreturn=current->num_descndnt_interior;
//}



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


//int enumerate_internal_branch(Node *current,int e_num_old){
/*! \brief enumerate the internal branches */
int Net::enumerate_internal_branch(Node *current, /*!< pointer to the node that is enumerated */
	int e_num_old) /*!< enumerator which is about to be updated \todo change e_num_old to int* type */ 
{
//needs modification, this is not correct.	
	
	if (debug_bool){
		cout<<"Net::enumerate_internal_branch start"<<endl;
	}
	int e_num_new;
	if (current->tip_bool){
		e_num_new=e_num_old;}
	else {
		if ((current->visited)==true){
			e_num_new=e_num_old+1;
			current->e_num2=e_num_new;
			}
		else{
			int e_num_dummy=e_num_old;
			for (int i_num_child=0;i_num_child<current->num_child;i_num_child++){
				e_num_dummy=enumerate_internal_branch(current->child[i_num_child],e_num_dummy);
			}
			current->visited=true;
			e_num_new=e_num_dummy+1;
			//if (current->e_num>0 && ){
				//current->e_num2=e_num_new;
			//}
			//else{
			current->e_num=e_num_new;
			//}
		}
	}
	if (debug_bool){
		cout<<"Net::enumerate_internal_branch end"<<endl;
	}
	return e_num_new;
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


/*! \brief free the meomory */
void Net::clear(){
	tax_name.clear();
	Net_nodes.clear();
};


bool Net::is_net_func(){
	//false stands for tree, true stands for net_work
	bool is_net_return=false;
	for (unsigned int i=0;i<Net_nodes.size();i++){
		if (Net_nodes[i].parent2){
			is_net_return=true;
			break;
		}
	}
	return is_net_return;	/*!< if its net or not*/
}

void Net::print_all_node(){
	//bool debug_switch=false;
	if ( is_net ){
		cout<<"           label  hybrid hyb_des non-tp parent1  abs_t brchln1 parent2 brchln2 #child #dsndnt #id rank   e_num   Clade "<<endl;
		//if (debug_switch){
			//string bug_statement="           label  hybrid hyb_des non-tip parent1 absolute_t brchlen1 parent2 brchlen2 #child #dsndnt #id rank   e_num   Clade ";
			////appending_debug_file(bug_statement);
		//}
		for (unsigned int i=0;i<Net_nodes.size();i++){
			for (unsigned int j=0;j<descndnt[i].size();j++){
				cout<<descndnt[i][j];
			}

			Net_nodes[i].print_net_Node();
			cout<<"  ";
			for (unsigned int j=0;j<descndnt2[i].size();j++){
				cout<<descndnt2[i][j];
			}
			cout<<endl;
		}	
	}
	else{
		cout<<"            label non-tp   parent        abs_t brchln #child #dsndnt #id rank e_num   Clade "<<endl;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			for (unsigned int j=0;j<descndnt[i].size();j++){
				cout<<descndnt[i][j];
			}
			Net_nodes[i].print_tree_Node();
						cout<<"  ";
			for (unsigned int j=0;j<descndnt2[i].size();j++){
				cout<<descndnt2[i][j];
			}
			cout<<endl;
		}	
	}		
}

/*! \brief Label interior node if the interior nodes of the tree string are not labeled */
string Net::label_interior_node(string in_str /*!< input newick form string */){
	vector <string> in_str_partition;
	int interior_node_counter=0;
	int sub_str_start_index=0;			
	size_t i=in_str.find(')');
	while ( i<in_str.size() ){
		interior_node_counter++;
		string current_string;
		size_t found_next_bracket=min(in_str.find(")",sub_str_start_index),in_str.size());
		current_string=in_str.substr(sub_str_start_index,found_next_bracket - sub_str_start_index +1);
		if (in_str[i+1]==';' || i==(in_str.size()-1)){
			current_string=current_string+"root";
			in_str_partition.push_back(current_string);
		}
		else{
			ostringstream interior_node_counter_str;
			interior_node_counter_str<<interior_node_counter;
			current_string=current_string+"Int_";
			current_string=current_string+interior_node_counter_str.str();
			in_str_partition.push_back(current_string);
			sub_str_start_index=i+1;
		}
		i=in_str.find(")",i+1);
	}
	string out_str;
	for (size_t i=0;i<in_str_partition.size();i++){
		//cout<<in_str_partition[i]<<endl;
		out_str=out_str+in_str_partition[i];
	}
	return out_str;
}

void appending_debug_file(string debug_file_input){
	ofstream debug_file;
	debug_file.open ("debug_file", ios::out | ios::app | ios::binary); 
	debug_file << debug_file_input << "\n";
	debug_file.close();
}

/*! \brief Add more information to log_file */
void appending_log_file(string log_file_input /*! Information added*/){
	ofstream log_file;
	log_file.open ("log_file", ios::out | ios::app | ios::binary); 
	log_file << log_file_input << "\n";
	log_file.close();
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

/*! \brief Remove interior nodes label of a string */
string remove_interior_label(string in_str/*!< input newick form string */){
	string out_str;
	out_str=in_str;
	//for (size_t char_i=0;char_i<out_str.size();){
		//if (out_str[char_i]==')' && (isalpha(out_str[char_i+1]) || isdigit(out_str[char_i+1]))){
		////if (out_str[char_i]==')' && isalpha(out_str[char_i+1]) ){
			////size_t char_j;
			////for (char_j=char_i+1;char_j<out_str.size();char_j++){
				////if (out_str[char_j]==':' || out_str[char_j]==',' || out_str[char_j]==';' || out_str[char_j]==')'){		
					////break;
				////}
			////}
			////cout<<out_str[char_j-1]<<out_str[char_j] <<"  "<<"working "<<char_j<<" not working "<<end_of_label_or_bl(out_str, char_i+1) +1<<endl;
			////out_str.erase(out_str.begin()+char_i+1,out_str.begin()+char_j);	
			//size_t char_j=end_of_label_or_bl(out_str, char_i+1);
			//out_str.erase(out_str.begin()+char_i+1,out_str.begin()+char_j+1);	
		//}
		//else{
			//char_i++;
		//}
	//}
	
	size_t found_bracket=out_str.find(')');
	while ( found_bracket<out_str.size() ){
		if (isalpha(out_str[found_bracket+1]) || isdigit(out_str[found_bracket+1])){
			size_t char_j=end_of_label_or_bl(out_str, found_bracket+1);
			//cout<<out_str<<endl;
			//cout<<out_str[char_j+1]<<endl;
			out_str.erase(out_str.begin()+found_bracket+1,out_str.begin()+char_j+1);
			//cout<<out_str<<endl;
		}
		found_bracket=out_str.find(")",found_bracket+1);
	}

	return out_str;
}


/*! \brief Remove branch length in a tree string, gives the tree topology */
string remove_brchlen(string in_str /*!< input newick form string */){
	string out_str = in_str;
	size_t found_col=out_str.find(':');
	while ( found_col<out_str.size() ){
		size_t char_j=end_of_label_or_bl(out_str,found_col)+1;
		//out_str.erase(out_str.begin()+found_col,out_str.begin()+char_j);
		out_str.erase(found_col,char_j-found_col);
		found_col=out_str.find(":",found_col+1);
	}
	return out_str;
}

/*! \brief Remove branch length in a tree string, gives the tree topology \todo combine with remove_brchlen, only keep one of them*/
string tree_topo(string in_str /*!< input newick form string */){
	return remove_brchlen(remove_interior_label(in_str));
}

string construct_adding_new_Net_str(Net in_Net){
	string out_str;
	out_str=in_Net.Net_nodes.back().node_content;
	out_str=out_str+in_Net.Net_nodes.back().label;
	if (in_Net.Net_nodes.back().brchlen1!=0){
		ostringstream brchlen_str;
		brchlen_str<<in_Net.Net_nodes.back().brchlen1;
		out_str=out_str+":"+brchlen_str.str();
	}
	out_str.push_back(';');
	return out_str;
}

/*! \brief Produce a tex file, which is used to draw the network 
 */
void plot_in_latex_file(const char* file_name /*! Name for the figure file */ , 
	Net net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/ 
	){
	ofstream latex_file;
	string file_name_no_dot(file_name);
	string file_name_with_dot(file_name);
	string dtex(".tex");
	
	size_t found=file_name_no_dot.find(dtex);
	if (found!=string::npos){
		file_name_no_dot=file_name_no_dot.substr(0,found);
	}
	else{
		file_name_with_dot=file_name_with_dot+dtex;
	}
	
	latex_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
	latex_file <<"\\documentclass[10pt]{article}\n";
	latex_file <<"\\usepackage{tikz,graphics,graphicx,lscape,fullpage,multicol,setspace}\n \\singlespacing\n \\begin{document}\n ";	
	latex_file<<"\\ifx\\du\\undefined\\newlength{\\du}\\fi\\setlength{\\du}{30\\unitlength}\n";
	latex_file <<"\\begin{center}\n";
	latex_file.close();
	plot_in_latex(file_name_with_dot.c_str(), net_dummy,plot_option);	
	//		plot_in_latex(file_name, net_str,plot_option);	
	latex_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
	latex_file <<"\\end{center}\n";
	latex_file <<"\\end{document}\n";
	latex_file.close();
	
	string command="pdflatex "+file_name_no_dot+".tex";
	int sys=system(command.c_str());

	string appending_log_str="Network figure generated in file: "+file_name_no_dot+".pdf";
	appending_log_file(appending_log_str);

}



/*! \brief Core function of drawing a network in .tex files. 
 */
void plot_in_latex(const char* file_name /*! Name for the figure file */ , 
	Net net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/
	){
	//Net net_dummy(net_str);
	ofstream latex_file;
	latex_file.open (file_name, ios::out | ios::app | ios::binary); 	
	latex_file <<"\\begin{tikzpicture}[thick]\n";
	valarray <int>  x_node=det_x_node (net_dummy);
	for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size();node_i++){
		string sp_node_label=net_dummy.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		if (net_dummy.Net_nodes[node_i].tip_bool){
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<net_dummy.Net_nodes[node_i].rank<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		//latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<y_node[node_i]<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
		else{
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<net_dummy.Net_nodes[node_i].rank<<"\\du) [circle,fill=orange,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
	}

	for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size()-1;node_i++){
		string sp_node_label=net_dummy.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		string sp_node_parent1_label=net_dummy.Net_nodes[node_i].parent1->label;
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!net_dummy.Net_nodes[node_i].tip_bool){
			if (plot_option==1){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].e_num <<"};\n";
			}
			else{
				if (plot_option==2){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen1 <<"};\n";
				}
				else{
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
				}	
			}
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;	
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==1){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].e_num2 <<"};\n";
				}
				else{
					if (plot_option==2){
						latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen2 <<"};\n";
					}
					else{
						latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<");\n";
					}		
				}
			}
		}
		else{
			if (plot_option==2){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen1 <<"};\n";
			}
			else{
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
			}	
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==2){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen2 <<"};\n";
				}
				else{
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<");\n";
				}		
			}
		}
	}
	latex_file <<"\\end{tikzpicture}\n\n";
	latex_file.close();	
}

/*! \brief When drawing network in .tex files, detemine the x coordinates of nodes
 */
valarray <int>  det_x_node (Net net_dummy){
	valarray <int>  x_node (net_dummy.Net_nodes.size());
	x_node[x_node.size()-1]=0;
	
	//valarray <int>  y_node (net_dummy.Net_nodes.size());
	//y_node[x_node.size()-1]=net_dummy.Net_nodes.back().e_num;
	
	for (int rank_i=net_dummy.Net_nodes.back().rank;rank_i>0;rank_i--){
		vector <int> x_node_dummy;
		vector <unsigned int> x_node_dummy_index;
		for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size();node_i++){
			if (net_dummy.Net_nodes[node_i].rank==rank_i){
				unsigned int n_child=net_dummy.Net_nodes[node_i].child.size();
				int parent_x=x_node[node_i];
				//int parent_y=y_node[node_i];
				int start_child_x=parent_x-floor(n_child/2);
				//int start_child_x=0;
				bool odd_num_child=false;
				if ((n_child % 2) == 1){
					odd_num_child=true;
				}
				if (odd_num_child){
					for (unsigned int child_i=0;child_i<n_child;child_i++){
						for (unsigned int node_j=0;node_j<net_dummy.Net_nodes.size();node_j++){
							if (net_dummy.Net_nodes[node_j].label==net_dummy.Net_nodes[node_i].child[child_i]->label){			
								//int child_y=y_node[node_j];
								if (start_child_x==parent_x){
									x_node[node_j]=parent_x;
									//y_node[node_j]=parent_y-1;
									start_child_x++;
								}
								else{
									//x_node[node_j]=start_child_x*(parent_y-child_y)+parent_x;
									x_node[node_j]=start_child_x;
									//y_node[node_j]=parent_y-1;
								}
								start_child_x++;
							}
						}
					}
				}
				else{
					for (unsigned int child_i=0;child_i<n_child;child_i++){
						for (unsigned int node_j=0;node_j<net_dummy.Net_nodes.size();node_j++){
							if (net_dummy.Net_nodes[node_j].label==net_dummy.Net_nodes[node_i].child[child_i]->label){
								if (start_child_x==parent_x){										
									start_child_x++;
								}
								//	x_node[node_j]=start_child_x*(parent_y-child_y)+parent_x;
								x_node[node_j]=start_child_x;
								//y_node[node_j]=parent_y-1;
								start_child_x++;
							}
						}
					}
					
				}
				x_node_dummy.push_back(x_node[node_i]);
				x_node_dummy_index.push_back(node_i);
			}
		}
		if (x_node_dummy.size()>1){
			bool need_to_shift=true;
			while (need_to_shift){
				for (unsigned int x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (unsigned int x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
						if (current_x_node_dummy==x_node_dummy[x_node_dummy_j]){
							if (x_node_dummy[x_node_dummy_j]>0){
								x_node_dummy[x_node_dummy_j]++;
								x_node[x_node_dummy_index[x_node_dummy_j]]++;
							}
							else{
								x_node_dummy[x_node_dummy_j]--;
								x_node[x_node_dummy_index[x_node_dummy_j]]--;
							}
						}
					}
				}
				need_to_shift=false;
				for (unsigned int x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (unsigned int x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
						if (current_x_node_dummy==x_node_dummy[x_node_dummy_j]){
							need_to_shift=true;
							break;
						}		
					}
					if (need_to_shift){
						break;
					}
				}
			}
		}
	}
	return x_node;
	
}

/*! \brief Produce a dot file, which is used to draw the network, and compile the dot file to a pdf file.
 */
void plot_in_dot(const char* file_name /*! Name for the figure file */,
	Net net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/){
	string file_name_no_dot(file_name);
	string file_name_with_dot(file_name);
	string ddot(".dot");
	
	size_t found=file_name_no_dot.find(ddot);
	if (found!=string::npos){
		file_name_no_dot=file_name_no_dot.substr(0,found);
	}
	else{
		file_name_with_dot=file_name_with_dot+ddot;
	}

	//Net net_dummy(net_str);
	ofstream dot_file;
	check_and_remove(file_name_with_dot.c_str());
	dot_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
			
	dot_file <<"graph G {\n rankdir=BT; ratio=compress;\n";//page="14,14"; determines the size of the ps output
	//valarray <int>  x_node=det_x_node (net_dummy);

	for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size()-1;node_i++){
		string sp_node_label=net_dummy.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		string sp_node_parent1_label=net_dummy.Net_nodes[node_i].parent1->label;
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!net_dummy.Net_nodes[node_i].tip_bool){		
			if (plot_option==1){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num <<"\"];\n";
			}
			else{
				if (plot_option==2){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen1 <<"\"];\n";	
				}
				else{
					dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
				}	
			}
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==1){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num2 <<"\"];\n";
				}
				else{
					if (plot_option==2){
						dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen2 <<"\"];\n";	
					}
					else{
						dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num2 <<"\"];\n";
					}	
				}	
				
			}
		}
		else{
			if (plot_option==2){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen1 <<"\"];\n";	
			}
			else{
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
			}	
			
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);					
				if (plot_option==2){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen2 <<"\"];\n";	
				}
				else{
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num2 <<"\"];\n";
				}
			}
		}
	}
	
	if (net_dummy.is_ultrametric){
		for (int rank_i=net_dummy.Net_nodes.back().rank;rank_i>0;rank_i--){
			dot_file<<"{ rank=same; ";
			vector <int> x_node_dummy;
			vector <unsigned int> x_node_dummy_index;
			for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size();node_i++){
				if (net_dummy.Net_nodes[node_i].rank==rank_i){
					string sp_node_label=net_dummy.Net_nodes[node_i].label;
					sp_node_label=rm_and_hash_sign(sp_node_label);
					dot_file<<sp_node_label<<" ";
				}	
			}
			dot_file<<"} ;\n";
		}
	}
	dot_file <<"}\n";
	dot_file.close();

	string command="dot -Tps "+file_name_no_dot+".dot -o "+file_name_no_dot+".ps";
	int sys=system(command.c_str());
	command="convert "+file_name_no_dot+".ps -resize 100\% "+file_name_no_dot+".pdf";
	//cout<<command<<endl;
	sys=system(command.c_str());

	string appending_log_str="Dot figure generated in file: "+file_name_no_dot+".pdf";
	appending_log_file(appending_log_str);
}


/*! \brief Compute factorial of a \return double a! */
double factorial (double a){
	if (a > 1){
		return (a * factorial (a-1));}
	else{
		return (1);}
}


/*! \brief Compute a permutations of n \return double */
double n_permu_a (double n, double a){
	if (a>1){
		return (n*n_permu_a(n-1,a-1));
	}
	else{
		if (a==1){
			return (n);
		}
		else{
			return (1);
		}
	}
}

/*! \brief Compute n choose k \return double */
double n_choose_k(double n, double k){
	if (k<(n/2)){
		return (n_choose_k(n,n-k));}
	else{
		return (n_permu_a(n,k)/factorial(k));}
}

/*! \brief Compute factorial of a \return int a! */
int factorial_int (int a){
	if (a > 1){
		return (a * factorial_int (a-1));}
	else{
		return (1);}
}

/*! \brief Compute a permutations of n \return int */
int n_permu_a_int (int n, int a){
	if (a>1){
		return (n*n_permu_a_int(n-1,a-1));}
	else{
		if (a==1){
			return (n);}
		else{
			return (1);}
	}
}

/*! \brief Compute n choose k \return int */
int n_choose_k_int(int n, int k){
	if (k<(n/2)){
		return (n_choose_k_int(n,n-k));}
	else{
		return (n_permu_a_int(n,k)/factorial_int(k));}
}

/*! \brief Remove the '&' and '#' signs from a string \return string */
string rm_and_hash_sign(string in_str){
	//string out_str=in_str;
	//if (int(out_str.find('#'))>0){
		//out_str=out_str.substr(0,out_str.find('#'));
	//}
	//if (int(out_str.find('&'))>0){
		//out_str=out_str.substr(0,out_str.find('&'));
	//}
	//return out_str;	
	return rm_hash_sign(rm_and_sign(in_str));
}


/*! \brief Remove '#' signs and the gamma parameter from a string \return string */
string rm_hash_sign(string in_str){
	while (int(in_str.find('#'))>0 && in_str.find('#')!=string::npos) {
		size_t i=end_of_label_or_bl(in_str, in_str.find('#'));
		in_str.erase(in_str.find('#'),i-in_str.find('#')+1);
		//in_str.erase(in_str.find('#'),1);
	}
	return in_str;
}

string rm_and_sign(string in_str){
	while (int(in_str.find('&'))>0 && in_str.find('&')!=string::npos) {
		in_str.erase(in_str.find('&'),1);
	}
	return in_str;
}

//string rm_and_sign(string old_str){
	//string new_str=old_str;
	//int new_str_len=new_str.size();
	//for (int i_str_len=0;i_str_len<new_str_len;){
		//if (isalpha(new_str[i_str_len])){
			////int str_start_index=i_str_len;
			//int j_str_len;
			//vector <string> label;
			//for (j_str_len=i_str_len;j_str_len<new_str_len;j_str_len++){
				//string label_dummy;
				//if (new_str[j_str_len] != '&'){
					//label_dummy=new_str[j_str_len];
					//label.push_back(label_dummy);
				//}

				//char stop=new_str[j_str_len+1];
				//if (stop==',' || stop==')' || stop==':' || stop==';'){
					//break;
				//}
			//}
			//unsigned int erase_len=j_str_len-i_str_len;
			//if (1<label.size() && label.size()<(erase_len+1)){
				//sort(label.begin(),label.end());
				//string insert_str;
				//for (unsigned int label_i=0;label_i<label.size();label_i++){
					//insert_str=insert_str+label[label_i];
				//}
				
				//new_str.erase(i_str_len,erase_len+1);
				////cout<<insert_str<<endl;
				//new_str.insert(i_str_len,insert_str);
				//i_str_len=i_str_len+insert_str.size();
			//}
			//else{
				//i_str_len=j_str_len+1;
			//}
		//}
		//else {i_str_len++;}
	//}
	//return new_str;	
//}

/*! \brief Terminate the program and print out the log file*/
int my_exit(){ // change my_exit() to return 0, and terminate program
	int sys=system("cat log_file");		
	//exit(1);
	return 0;
}


/*! \brief Set plot options */
/*! \return int '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/ 
int set_plot_option(bool plot_label,bool plot_branch){
	int plot_option=0;
	if (plot_label){
		plot_option=1;
		appending_log_file("Internal branches are labelled by post-order tree traversal.");
	}
	else{
		if (plot_branch){
			plot_option=2;
			appending_log_file("Branch lengths are labelled.");
		}
	}
	
	//int plot_option;
	//if (plot_label){
		//plot_option=1;
		//appending_log_file("Internal branches are labelled by post-order tree traversal.");
	//}
	//else if (plot_branch){
			//plot_option=2;
			//appending_log_file("Branch lengths are labelled.");
		//}
		
	//else{
		//plot_option=0;
	//}
	return plot_option;
}

/*! Check and remove files*/
void check_and_remove(const char* file_name){
	ifstream my_file(file_name);
	if (my_file.good())
	{
	  remove(file_name);
	}
}

size_t end_of_label_or_bl(string in_str, size_t i){
	size_t j ;
	for (j=i;j<in_str.size();j++){
		char stop=in_str[j+1];
		if (stop==',' || stop==')' || stop==':' || stop==';'){
			break;
		}
	}
	return j;
}

string extract_label(string in_str, size_t i){
	size_t j=end_of_label_or_bl(in_str, i);
	//cout<<"i="<<i<<", j="<<j<<endl;
	return in_str.substr(i,j+1-i);
}


size_t hybrid_hash_index(string in_str){
	return in_str.find('#');
}

string extract_hybrid_label(string in_str){
	size_t hash_index=hybrid_hash_index(in_str);
	return in_str.substr(0,hash_index);
}

string extract_hybrid_para_str(string in_str){
	size_t hash_index=hybrid_hash_index(in_str);
	return in_str.substr(hash_index+1);//,in_str.size()-1);
}

double extract_hybrid_para(string in_str){
	double para;
	istringstream para_istr(extract_hybrid_para_str(in_str));
	para_istr>>para;
	return para;
}


string read_input_line(char inchar[]){
	ifstream in_file;
	string out_str;
	in_file.open(inchar);
	if (in_file.good()){
		getline (in_file,out_str);}
	else{
		string dummy_str(inchar);
		if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos){
		out_str=dummy_str;
		}else{
			cout<<"Error: check input '"<<inchar<<"'"<<endl;
		}
	}
	in_file.close();
			
return 	out_str;
}

vector <string> read_input_lines(char inchar[]){
	vector <string> out_vec;
	ifstream in_file;
	in_file.open(inchar);
	string out_str;
	if (in_file.good()){
		getline (in_file,out_str);
		while (out_str.size()>0){   
			out_vec.push_back(out_str);
			getline (in_file,out_str);
		}
	}	
	else{
		string dummy_str(inchar);
		if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos){
			out_str=dummy_str;
			out_vec.push_back(out_str);
		}else{
			cout<<"Error: check input '"<<inchar<<"'"<<endl;
		}
	}
	in_file.close();
	return out_vec;	
}

bool is_num(char inchar[]){
	bool is_num_return=true;
	string in_str(inchar);
	for (size_t i=0;i<in_str.size();i++){
		if (isalpha(in_str[i]) && in_str[i]!='e'){
			is_num_return=false;
			break;
		}
	}
	return is_num_return;
}

string read_input_para(char inchar[],string in_str){
	
	string out_str;
	if (is_num(inchar)){
		istringstream para_istrm(inchar);
		double para;
		para_istrm>>para;
		out_str=write_para_into_tree(in_str, para);
	}
	else{
		out_str=read_input_line(inchar);
	}
	return out_str;
}
		
		
/*! \brief Write a fixed parameter into a externed newick formatted network string*/
string write_para_into_tree(string in_str /*! Externed newick formatted network string*/, 
double para /*! Coalescent parameter or fixed population sizes */){
	if (in_str.size()==0){
		cout<<"Define input tree (network) first!"<<endl;
		exit(1);;
	}
	Net para_Net(in_str);
	vector <Node*> para_Net_node_ptr;
	if (debug_bool){
		cout<<"write_para_into_tree flag"<<endl;
	}
	for (unsigned int node_i=0;node_i<para_Net.Net_nodes.size();node_i++){
		Node* new_node_ptr=NULL;
        para_Net_node_ptr.push_back(new_node_ptr);
        para_Net_node_ptr[node_i]=&para_Net.Net_nodes[node_i];
		para_Net_node_ptr[node_i]->brchlen1=para;
		if (para_Net.Net_nodes[node_i].hybrid){
			para_Net_node_ptr[node_i]->brchlen2=para;
		}
	}
	if (debug_bool){
		cout<<"write_para_into_tree flag"<<para<<endl;
	}
	para_Net_node_ptr.back()->brchlen1=para;
	rewrite_node_content(para_Net_node_ptr);
	string para_string=construct_adding_new_Net_str(para_Net);
	return para_string;	
}
		
