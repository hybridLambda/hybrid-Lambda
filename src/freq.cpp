/*
 * hybrid_sim is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012 Sha (Joe) Zhu
 * 
 * This file is part of hybrid_sim.
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
/*! \file freq.cpp
 * \brief Count frequencies of tree topologies */


#include"freq.hpp"

/*! \brief For given tree strings, differentiate topologies, and count frequencies for each topology
 * 
 */
topo_freq::topo_freq(vector <string> gt_strings){
	Net gt_dummy(gt_strings[0]);
	gene_topo.push_back(tree_topo(gt_strings[0]));
	gene_freq.push_back(1);
	for (unsigned int gt_string_i=1;gt_string_i<gt_strings.size();gt_string_i++){
		bool new_topo=true;
		for (unsigned int topo_i=0;topo_i<gene_topo.size();topo_i++){
			if (same_topo(gt_strings[gt_string_i], gene_topo[topo_i])){
				gene_freq[topo_i]=gene_freq[topo_i]+1;
				new_topo=false;
				break;
			}	
		}
		if (new_topo){
			gene_topo.push_back(tree_topo(gt_strings[gt_string_i]));
			gene_freq.push_back(1);
		}
	}	
}


/*! \brief determine two tree strings have the same topology or not
 */
bool same_topo(string gt_string1,string gt_string2){
	bool same_topo_return=false;
	Net gt1(gt_string1);
	Net gt2(gt_string2);

	if (gt1.tax_name.size()!=gt2.tax_name.size() || gt1.descndnt2.size()!=gt2.descndnt2.size()){
		same_topo_return=false;
	}
	else{
		bool same_tax_names=true;
		for (unsigned int tax_i=0;tax_i<gt1.tax_name.size();tax_i++){
			if (gt1.tax_name[tax_i]!=gt2.tax_name[tax_i]){
				//cout<<gt1.tax_name[tax_i]<<" "<<gt2.tax_name[tax_i]<<endl;
				same_tax_names=false;
				break;
			}
		}
		
		if (same_tax_names){
			//valarray <int> gt1_in_gt2(gt1.descndnt.size(),0);
			//valarray <int> gt2_in_gt1(gt2.descndnt.size(),0);
			vector <int> gt1_in_gt2(gt1.descndnt2.size(),0);
			vector <int> gt2_in_gt1(gt2.descndnt2.size(),0);
			same_topo_return=true;
			for (unsigned int node1_i=0;node1_i<gt1.descndnt2.size();node1_i++){
				for (unsigned int node2_i=0;node2_i<gt2.descndnt2.size();node2_i++){
					valarray <bool> comp=(gt1.descndnt2[node1_i]==gt2.descndnt2[node2_i]);
					if (comp.min()==1){
						//cout<<gt1.Net_nodes[node1_i].node_content<<endl;
						//cout<<gt2.Net_nodes[node2_i].node_content<<endl;
						gt1_in_gt2[node1_i]=1;
						gt2_in_gt1[node2_i]=1;
						//break;
					}
					
				}
			}
			//for (unsigned int node2_i=0;node2_i<gt2.descndnt2.size();node2_i++){
				//for (unsigned int node1_i=0;node1_i<gt1.descndnt2.size();node1_i++){
					//valarray <bool> comp=(gt1.descndnt2[node1_i]==gt2.descndnt2[node2_i]);
					//if (comp.min()==1){
						////cout<<gt1.Net_nodes[node1_i].node_content<<endl;
						////cout<<gt2.Net_nodes[node2_i].node_content<<endl;
						//gt1_in_gt2[node1_i]=1;
						//gt2_in_gt1[node2_i]=1;
						////break;
					//}
					
				//}
			//}
			for (unsigned int node1_i=0;node1_i<gt1.descndnt2.size();node1_i++){
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
	//cout<<gt_string1<<"  "<<gt_string2<<endl;
	//if (same_topo_return){cout<<"same"<<endl;}
	
	return same_topo_return;

}


/*! \brief Determine the RF distance of two trees
 * \todo UNTESTED!!!!!
 */
int RF_dist(string gt_string1,string gt_string2){
	int RF_dist_return=0;
	Net gt1(gt_string1);
	Net gt2(gt_string2);
	valarray <int> gt1_in_gt2(gt1.descndnt.size(),0);
	valarray <int> gt2_in_gt1(gt2.descndnt.size(),0);
	for (unsigned int node1_i=0;node1_i<gt1.descndnt2.size();node1_i++){
		for (unsigned int node2_i=0;node2_i<gt2.descndnt2.size();node2_i++){
			valarray <bool> comp=(gt1.descndnt2[node1_i]==gt2.descndnt2[node2_i]);
			if (comp.min()==1){
				//cout<<gt1.Net_nodes[node1_i].node_content<<endl;
				//cout<<gt2.Net_nodes[node2_i].node_content<<endl;
				gt1_in_gt2[node1_i]=1;
				gt2_in_gt1[node2_i]=1;
				//break;
			}
		}
	}
	RF_dist_return=gt1_in_gt2.size()-gt1_in_gt2.sum()+gt2_in_gt1.size()-gt2_in_gt1.sum();
	return RF_dist_return;

}


/*! \brief Compute gene tree frequencies */
void compute_gt_frequencies(vector <string> gt_tree_str_s, string freq_file_name){
	remove(freq_file_name.c_str());
	topo_freq my_topo_freq(gt_tree_str_s);
	ofstream freq_out_file;
	freq_out_file.open (freq_file_name.c_str(), ios::out | ios::app | ios::binary); 
	int total_num=0;
	for (unsigned int topo_i=0;topo_i<my_topo_freq.gene_topo.size();topo_i++){
		freq_out_file<<topo_i+1<<" "<<my_topo_freq.gene_topo[topo_i]<<"  "<<my_topo_freq.gene_freq[topo_i]<<endl;
		total_num=total_num+my_topo_freq.gene_freq[topo_i];
	}
	cout<<total_num<<endl;
	freq_out_file.close();
	string appending_log_str="Gene trees frequency analyzed in file: "+freq_file_name;
	appending_log_file(appending_log_str);
}




