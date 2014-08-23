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

/*! \file utility.cpp
 */

#include"utility.hpp"

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
		throw std::invalid_argument(in_str + "Parenthesis not balanced!" );
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


/*! \brief Remove interior nodes label of a string */
string remove_interior_label(string in_str/*!< input newick form string */){
	string out_str;
	out_str=in_str;
	
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
	
void readNextStringto( string &readto , int& argc_i, int argc_, char * const* argv_ ){
    argc_i++;
    if (argc_i >= argc_) throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv_[argc_i-1]); 
    readto = std::string(argv_[argc_i]);
    if ( readto[0] == '-' ) throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv_[argc_i-1]);
}
