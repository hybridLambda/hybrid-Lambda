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
//seg-site.cpp 

#include"seg-site.hpp"
#include"net.hpp"


seg::param::param(){
	seg_dir_name="seg-sites";
}

seg::param::param(int argc, char *argv[]){	
	seg_dir_name="seg-sites";
	for (int argc_i=1;argc_i < argc;argc_i++){
		std::string argv_i(argv[argc_i]);
		if (argv_i=="-segD"){
			seg_dir_name=argv[argc_i+1];		
			break;
		}	
	}
	check_and_remove(seg_dir_name.c_str());			
}


/*! \brief remove old segregating sites data, and generate new ones
 * \todo user defined direcotry name*/
void seg::param::create_site_data_dir(vector <string> mt_tree_str_s){
	string rm_commond="rm -rf "+ seg_dir_name;
	int sys=system(rm_commond.c_str());
	string mkdir_commond="mkdir "+seg_dir_name;
	 sys=system(mkdir_commond.c_str());
	for (size_t num_mt_tree_i=0;num_mt_tree_i<mt_tree_str_s.size();num_mt_tree_i++){
		create_new_site_data(mt_tree_str_s[num_mt_tree_i],num_mt_tree_i+1);
	}
}

/*! \brief Generate segrateing site data 
 * \todo user defined direcotry name */
void seg::param::create_new_site_data(string gt_string_mut_num,int site_i){
	Net mt_tree(gt_string_mut_num);
	ofstream site_data_file;
	ostringstream site_i_str;
	site_i_str<<site_i;
	string sitefile_name=seg_dir_name+"/site"+site_i_str.str();
	site_data_file.open (sitefile_name.c_str()); 
	//cout<<mt_tree.tip_name.size()<<endl;
		//cout<<mt_tree.tax_name.size()<<endl;
	
	int total_mut=0;
	for (unsigned int node_i=0;node_i<mt_tree.Net_nodes.size();node_i++){
		total_mut=total_mut+mt_tree.Net_nodes[node_i].brchlen1;
	}
	//site_data_file<<gt_string_mut_num<<endl;
	//site_data_file<<"Total "<<total_mut<<" mutations"<<endl<<endl;
	for (unsigned int tip_i=0;tip_i<mt_tree.tip_name.size();tip_i++){
		site_data_file<<mt_tree.tip_name[tip_i]<<" ";
		//cout<<mt_tree.tip_name[tip_i]<<" ";
		for (unsigned int node_i=0;node_i<mt_tree.Net_nodes.size();node_i++){
			if (mt_tree.Net_nodes[node_i].brchlen1>0){
				for (int num_repeat=0;num_repeat<mt_tree.Net_nodes[node_i].brchlen1;num_repeat++ ){				
					//site_data_file<<mt_tree.descndnt2[node_i][tip_i] <<" ";
					site_data_file<<mt_tree.descndnt2[node_i][tip_i] ;
				}
			}
		}
		site_data_file<<"\n";
	}
	site_data_file.close();	
}



