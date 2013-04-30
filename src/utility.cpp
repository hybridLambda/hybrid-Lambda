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

/*! \file utility.cpp
 *  \brief Core function of converting a Newick (extended Newick) format string into a species tree (network), and simple string manipulation for tree strings
 */



#include"utility.hpp"
//extern bool debug_bool;
//extern bool log_file_bool;



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




void appending_debug_file(string debug_file_input){
	ofstream debug_file;
	debug_file.open ("debug_file", ios::out | ios::app | ios::binary); 
	debug_file << debug_file_input << "\n";
	debug_file.close();
}

/*! \brief Add more information to log_file */
void appending_log_file(std::string log_file_NAME,std::string log_file_input /*! Information added*/){
	std::ofstream log_file;
	log_file.open (log_file_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
	log_file << log_file_input << "\n";
	log_file.close();
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
			//cout<<"Error: check input '"<<inchar<<"'"<<endl;
			string error_msg(inchar);
			
			throw std::invalid_argument("Invalid input file. " +error_msg);
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
			//cout<<"Error: check input '"<<inchar<<"'"<<endl;
			string error_msg(inchar);
			throw std::invalid_argument("Invalid input file. " +error_msg);
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
		
	
		
