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

#include "freq.hpp"

Frequency::Frequency( int argc, char * const* argv ):
    argc_(argc), argv_(argv){
	this->freq_out_filename = "freq_out";
	for (  argc_i = 1; argc_i < argc; argc_i++ ){
		std::string argv_i( argv[argc_i] );
		if ( argv_i == "-freq_file" || argv_i=="-fF" ){
            readNextStringto( this->freq_out_filename , this->argc_i, this->argc_,  this->argv_ ); 
			break;
		}
        argc_i++;	
	}

    ifstream tmp_file( this->freq_out_filename.c_str() );
	if ( tmp_file.good() ) 	{  remove(freq_out_filename.c_str()); 	}
}


/*! \brief For given tree strings, differentiate topologies, and count frequencies for each topology
 * 
 */  
/*! \brief Compute gene tree frequencies */
void Frequency::compute_gt_frequencies( vector <string> &gt_tree_str_s ){
    this->gt_tree_str_tmp = gt_tree_str_s;
	this->compute_gt_frequencies_core();
	freq_out_file.open (this->freq_out_filename.c_str(), ios::out | ios::app | ios::binary); 
	int total_num = 0;
	for (size_t topo_i = 0; topo_i < this->gene_topo.size(); topo_i++){
		freq_out_file << topo_i+1 << " " << this->gene_topo[topo_i] << "  " << this->gene_freq[topo_i] << endl;
		total_num = total_num + this->gene_freq[topo_i];
	}
	dout << total_num << endl;
	freq_out_file.close();
    clog << "Frequency file is saved at: " << this->freq_out_filename << "\n";
}

void Frequency::compute_gt_frequencies_core(){
    Net checking_topo( this->gt_tree_str_tmp[0] );
	this->gene_topo.push_back( tree_topo( this->gt_tree_str_tmp[0]) );
	this->gene_freq.push_back(1);
	for (size_t i = 1; i < this->gt_tree_str_tmp.size(); i++ ){
		bool new_topo=true;
		for (size_t topo_i = 0; topo_i < gene_topo.size(); topo_i++ ){
			if (same_topo(this->gt_tree_str_tmp[i], gene_topo[topo_i])){
				gene_freq[topo_i]=gene_freq[topo_i]+1;
				new_topo = false;
				break;
			}	
		}
		if ( new_topo ){
			gene_topo.push_back( tree_topo(this->gt_tree_str_tmp[i]) );
			gene_freq.push_back(1);
		}
	}	
}


/*! \brief determine two tree strings have the same topology or not
 */
bool Frequency::same_topo( string gt_string1, string gt_string2 ){
	bool same_topo_return=false;
	Net gt1(gt_string1);
	Net gt2(gt_string2);

	if (gt1.tax_name.size()!=gt2.tax_name.size() || gt1.samples_below.size()!=gt2.samples_below.size()){
		same_topo_return=false;
	}
	else{
		bool same_tax_names=true;
		for (size_t tax_i=0;tax_i<gt1.tax_name.size();tax_i++){
			if (gt1.tax_name[tax_i]!=gt2.tax_name[tax_i]){
				dout<<gt1.tax_name[tax_i]<<" "<<gt2.tax_name[tax_i]<<endl;
				same_tax_names=false;
				break;
			}
		}
		
		if (same_tax_names){
			vector <int> gt1_in_gt2(gt1.samples_below.size(),0);
			vector <int> gt2_in_gt1(gt2.samples_below.size(),0);
			same_topo_return=true;
			for (size_t node1_i=0;node1_i<gt1.samples_below.size();node1_i++){
				for (size_t node2_i=0;node2_i<gt2.samples_below.size();node2_i++){
					valarray <bool> comp=(gt1.samples_below[node1_i]==gt2.samples_below[node2_i]);
					if (comp.min()==1){
						gt1_in_gt2[node1_i]=1;
						gt2_in_gt1[node2_i]=1;
					}					
				}
			}
			for (size_t node1_i=0;node1_i<gt1.samples_below.size();node1_i++){
				if (gt1_in_gt2[node1_i]==0 || gt2_in_gt1[node1_i]==0){
					same_topo_return=false;
					break;// this is newly added. test????
				}
			}
		}
		else{
			same_topo_return=false;
		}
	}
	if (same_topo_return){dout<<"same"<<endl;}	
	return same_topo_return;
}

/*! \brief Remove branch length in a tree string, gives the tree topology */
string Frequency::remove_brchlen(string in_str /*!< input newick form string */){
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
string Frequency::tree_topo(string in_str /*!< input newick form string */){
	return remove_brchlen(remove_interior_label(in_str));
}

