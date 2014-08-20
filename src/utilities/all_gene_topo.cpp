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

/*! \file all_gene_topo.cpp
 * 	\brief  Enumerate all possible tree topologies for given number of taxon */

#include"all_gene_topo.hpp"


vector < string > find_current_taxa_name( string in_str ){
	vector < string > current_taxa_name;
	for (size_t i=1;i<in_str.size();){
		if (isalpha(in_str[i]) || isdigit(in_str[i])){
			string label=extract_label(in_str, i);
			current_taxa_name.push_back(label);
			i=label.size()+i;
		}
		else {
			i++;
		}
	}
	return current_taxa_name;
}


string add_new_taxa_at_tip(string in_str,size_t i, size_t tax_i, vector <string> taxa_name,string label){
	string out_str=in_str;			
	string new_add_in="("+label+","+taxa_name[tax_i]+")";				
	out_str.replace(i,label.size(),new_add_in);
	return out_str;
}


string add_new_taxa_at_int(string in_str, size_t i, size_t tax_i, vector <string> taxa_name){
	string out_str=in_str;
	size_t rev_dummy_i=Parenthesis_balance_index_backwards(in_str,i);
	size_t substr_len=i-rev_dummy_i+1;
	string new_add_in=in_str.substr(rev_dummy_i,substr_len);
	new_add_in="("+new_add_in+","+taxa_name[tax_i]+")";
	out_str.replace(rev_dummy_i,substr_len,new_add_in);
	return out_str;
	
}


/*! \fn vector <string> generate_new_topo_list(vector <string> old_topo_list, size_t tax_num)
* \brief Enumerate new tree topologies for k+1 taxon given tree topologies for k taxon 
* \return the list of \a n taxon topologies 
*/
vector <string> generate_new_topo_list(
vector <string> old_topo_list /*! topology list which is about to be updated */,
vector <string> taxa_name /*! list of taxa names */)
{
	//size_t tax_num=taxa_name.size();
	vector <string> new_topo_list;
	vector <string> out_topo_list;
	vector < string > current_taxa_name=find_current_taxa_name(old_topo_list[0]);
			
	if (current_taxa_name.size()<taxa_name.size()){
		for (size_t i=0;i<old_topo_list.size();i++){
			string old_topo_list_dummy=old_topo_list[i];
			for (size_t i_str_len=1;i_str_len<old_topo_list_dummy.size();){				
				if ( isalpha(old_topo_list_dummy[i_str_len]) || isdigit(old_topo_list_dummy[i_str_len])){
					 
					string label=extract_label(old_topo_list_dummy, i_str_len);
					//string new_topo_list_dummy=old_topo_list_dummy;			
					//string new_add_in="("+label+","+taxa_name[current_taxa_name.size()]+")";				
					//new_topo_list_dummy.replace(i_str_len,label.size(),new_add_in);
					
					string new_topo_list_dummy=add_new_taxa_at_tip(old_topo_list_dummy,i_str_len, current_taxa_name.size(), taxa_name,label);
					new_topo_list.push_back(new_topo_list_dummy);
					i_str_len=label.size()+i_str_len;
				}
				else{
					if (old_topo_list_dummy[i_str_len]==')'){
						//size_t rev_dummy_i=Parenthesis_balance_index_backwards(old_topo_list_dummy,i_str_len);
						//size_t substr_len=i_str_len-rev_dummy_i+1;
						//string new_topo_list_dummy=old_topo_list_dummy;					
						//string new_add_in=old_topo_list_dummy.substr(rev_dummy_i,substr_len);
						//new_add_in="("+new_add_in+","+taxa_name[current_taxa_name.size()]+")";
						//new_topo_list_dummy.replace(rev_dummy_i,substr_len,new_add_in);
						
						string new_topo_list_dummy=add_new_taxa_at_int(old_topo_list_dummy, i_str_len, current_taxa_name.size(), taxa_name);
						new_topo_list.push_back(new_topo_list_dummy);					
					}
					i_str_len++;
				}
			}
		}
	}
	
	if ((current_taxa_name.size()+1)<taxa_name.size()){
		out_topo_list=generate_new_topo_list(new_topo_list, taxa_name);
	}
	else{
		out_topo_list=new_topo_list;
	}

	return out_topo_list;
}


/*! \fn vector <string> all_n_tax_gene_tree(size_t tax_num)
 * \brief Enumerate all possible topologies for n taxon 
 * \return the list of \a n taxon topologies*/
vector <string> all_n_tax_gene_tree(vector <string> taxa_name/*! list of taxa names */){
	vector < string > old_topo_list;
	vector < string > new_topo_list;
	string first_topo="("+taxa_name[0]+","+taxa_name[1]+")";
	old_topo_list.push_back(first_topo);
	new_topo_list=generate_new_topo_list(old_topo_list, taxa_name);
	return new_topo_list;
}

