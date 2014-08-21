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
/*! \file freq.cpp
 * \brief Count frequencies of tree topologies */

#include"freq.hpp"

/*! \brief Determine the RF distance of two trees
 * \todo UNTESTED!!!!!
 */
int RF_dist(string gt_string1,string gt_string2){
	int RF_dist_return=0;
	Net gt1(gt_string1);
	Net gt2(gt_string2);
	valarray <int> gt1_in_gt2(gt1.descndnt.size(),0);
	valarray <int> gt2_in_gt1(gt2.descndnt.size(),0);
	for (size_t node1_i=0;node1_i<gt1.descndnt2.size();node1_i++){
		for (size_t node2_i=0;node2_i<gt2.descndnt2.size();node2_i++){
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
