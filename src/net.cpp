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

#include"net.hpp"

/*! \brief Construct Net object from a (extended) Newick string */
Net::Net(string old_string /*! input (extended) newick form string */){
	if (old_string.size()>0){

		checking_Parenthesis(old_string);
		net_str=checking_labeled(old_string);
	
		//dout<<"Net::Net flag1"<<endl;

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
			
		
	//		dout<<"Net::Net flag2"<<endl;
		
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
		
			//dout<<"Net::Net flag3"<<endl;
		
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
		
			//dout<<"Net::Net flag4"<<endl;
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
			//dout<<Net_nodes_ptr[i]->name<<"  "<<Net_nodes_ptr[i]->label<<"  "<<Net_nodes_ptr[i]->tip_bool << " ";
			for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				//dout<<Net_nodes_ptr[i]->label<<"   "<<tax_name[tax_name_i]<<"  ";
				if (find_descndnt(Net_nodes_ptr[i],tax_name[tax_name_i])){
					descndnt[i][tax_name_i]=1;
				}
				//dout<<descndnt[i][tax_name_i];
			}
			//dout<<endl;
			for (unsigned int tip_name_i=0;tip_name_i<tip_name.size();tip_name_i++){
				if (find_descndnt2(Net_nodes_ptr[i],tip_name[tip_name_i])){
					descndnt2[i][tip_name_i]=1;
				}
			}
			Net_nodes[i].num_descndnt=descndnt[i].sum();
		}
			//dout<<"Net::Net flag5"<<endl;
		
		//num_descndnt_interior(Net_nodes_ptr.back()); //replace by following
		for (unsigned int i=0;i<Net_nodes_ptr.size();i++){
			for (unsigned int j=0;j<Net_nodes_ptr.size();j++){
				if (i!=j){
					valarray <int> descndnt_diff=(descndnt[i]-descndnt[j]);
					if (descndnt_diff.min() >= 0 && Net_nodes[i].rank > Net_nodes[j].rank && Net_nodes[j].rank>=2){
						Net_nodes_ptr[i]->num_descndnt_interior=Net_nodes_ptr[i]->num_descndnt_interior+1;
						Net_nodes_ptr[i]->descndnt_interior_node.push_back(Net_nodes_ptr[j]);
					}
					
				}
			
			}
			//dout<<"checking #interior_des "<<Net_nodes_ptr[i]->num_descndnt_interior<<" "<<Net_nodes_ptr[i]->descndnt_interior_node.size()<<endl;
		}
			
		int e_num_old=0;
		enumerate_internal_branch(Net_nodes_ptr.back(),e_num_old);
		
			//dout<<"Net::Net flag6.5"<<endl;
		
		//dout<<descndnt.size()<<endl;
		//dout<<net_str<<endl;
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
	//dout<<"Net constructed"<<endl;
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



//int enumerate_internal_branch(Node *current,int e_num_old){
/*! \brief enumerate the internal branches */
int Net::enumerate_internal_branch(Node *current, /*!< pointer to the node that is enumerated */
	int e_num_old) /*!< enumerator which is about to be updated \todo change e_num_old to int* type */ 
{
	
		//dout<<"Net::enumerate_internal_branch start"<<endl;
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
		//dout<<"Net::enumerate_internal_branch end"<<endl;
	return e_num_new;
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
	if ( is_net ){
		cout<<"           label  hybrid hyb_des non-tp parent1  abs_t brchln1 parent2 brchln2 #child #dsndnt #id rank   e_num   Clade "<<endl;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			for (unsigned int j=0;j<descndnt[i].size();j++){
				cout<<setw(3)<<descndnt[i][j];
			}

			Net_nodes[i].print_net_Node();
			cout<<"  ";
			for (unsigned int j=0;j<descndnt2[i].size();j++){
				cout<<setw(3)<<descndnt2[i][j];
			}
			cout<<endl;
		}	
	}
	else{
		cout<<"            label non-tp   parent        abs_t brchln #child #dsndnt #id rank e_num   Clade "<<endl;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			for (unsigned int j=0;j<descndnt[i].size();j++){
				cout<<setw(3)<<descndnt[i][j];
			}
			Net_nodes[i].print_tree_Node();
						cout<<"  ";
			for (unsigned int j=0;j<descndnt2[i].size();j++){
				cout<<setw(3)<<descndnt2[i][j];
			}
			cout<<endl;
		}	
	}		
}

void Net::print_all_node_dout(){
	if ( is_net ){
		dout<<"           label  hybrid hyb_des non-tp parent1  abs_t brchln1 parent2 brchln2 #child #dsndnt #id rank   e_num   Clade "<<endl;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			for (unsigned int j=0;j<descndnt[i].size();j++){
				dout<<setw(3)<<descndnt[i][j];
			}

			//Net_nodes[i].print_net_Node();
			Net_nodes[i].print_net_Node_dout();
			dout<<"  ";
			for (unsigned int j=0;j<descndnt2[i].size();j++){
				dout<<descndnt2[i][j];
			}
			dout<<endl;
		}	
	}
	else{
		dout<<"            label non-tp   parent        abs_t brchln #child #dsndnt #id rank e_num   Clade "<<endl;
		for (unsigned int i=0;i<Net_nodes.size();i++){
			for (unsigned int j=0;j<descndnt[i].size();j++){
				dout<<setw(3)<<descndnt[i][j];
			}
			Net_nodes[i].print_tree_Node_dout();
						dout<<"  ";
			for (unsigned int j=0;j<descndnt2[i].size();j++){
				dout<<descndnt2[i][j];
			}
			dout<<endl;
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
		//dout<<in_str_partition[i]<<endl;
		out_str=out_str+in_str_partition[i];
	}
	return out_str;
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


	
/*! \brief Write a fixed parameter into a externed newick formatted network string*/
string write_para_into_tree(string in_str /*! Externed newick formatted network string*/, 
double para /*! Coalescent parameter or fixed population sizes */){
	if (in_str.size()==0){
		throw std::invalid_argument("Please define the input tree (network).");
	}
	Net para_Net(in_str);
	vector <Node*> para_Net_node_ptr;
	for (unsigned int node_i=0;node_i<para_Net.Net_nodes.size();node_i++){
		Node* new_node_ptr=NULL;
        para_Net_node_ptr.push_back(new_node_ptr);
        para_Net_node_ptr[node_i]=&para_Net.Net_nodes[node_i];
		para_Net_node_ptr[node_i]->brchlen1=para;
		if (para_Net.Net_nodes[node_i].hybrid){
			para_Net_node_ptr[node_i]->brchlen2=para;
		}
	}
	para_Net_node_ptr.back()->brchlen1=para;
	rewrite_node_content(para_Net_node_ptr);
	string para_string=construct_adding_new_Net_str(para_Net);
	return para_string;	
}


size_t Net::first_coal_rank(){
    size_t min_rank = (size_t)Net_nodes.back().rank;
    for (size_t i = 0 ; i < Net_nodes.size(); i++){
        if ( !Net_nodes[i].tip_bool ){
            min_rank = ( Net_nodes[i].rank < min_rank ) ?  Net_nodes[i].rank : min_rank ;
            }
        }
    return min_rank;
    }

size_t Net::first_coal_index (){    
    size_t min_rank = this->first_coal_rank();
    size_t dummy_index = this->Net_nodes.size()-1;
    double min_coal_time = this->Net_nodes[dummy_index].absolute_time;
    //cout<<"min_rank = "<<min_rank<<endl;
    for (size_t i = 0 ; i < Net_nodes.size(); i++){
        if ( this->Net_nodes[i].rank == min_rank &&  this->Net_nodes[i].absolute_time < min_coal_time ){
            dummy_index = i;
            min_coal_time = this->Net_nodes[dummy_index].absolute_time;
            }
        
        }
        //cout << "min_coal_time = " << min_coal_time  <<endl;
    return dummy_index;
    }
