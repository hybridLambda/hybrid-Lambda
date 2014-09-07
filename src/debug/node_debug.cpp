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

#include "../node.hpp"
void Node::print_net_Node_dout(){
	dout<<setw(12)<<label;
	dout<<setw(6)<<hybrid;
	dout<<setw(8)<<descndnt_of_hybrid;
	dout<<setw(5)<<tip_bool;
    dout<<setw (11)<< (parent1) ? (parent1->label) : "" ;

	dout<<setw (8)<<height;
	dout<<setw (8)<<brchlen1;
    dout<<setw (11)<< (parent2) ? (parent2->label) : "" ;

	dout<<setw (8)<<brchlen2;
	dout<<setw (7)<<num_child;
	dout<<setw (8)<<num_descndnt;
	dout<<setw(4)<<num_descndnt_interior;
//	dout<<setw (7)<<current.e_num;
//	dout<<setw (3)<<current.e_num2;
	dout<<setw (6)<<rank<<"   ";
	//for (size_t i=0;i<descndnt.size();i++){
		//dout<<setw (1)<<descndnt[i];
	//}
	dout<<setw(2)<<e_num;
	dout<<setw(3)<<e_num2;
	dout<<"    "<<clade;
	//dout<<endl;
}


void Node::print_tree_Node_dout(){
	dout<<setw (12)<<label;
	dout<<setw(5)<<tip_bool;
    dout<<setw (11)<< (parent1) ? (parent1->label) : "" ;
	dout<<setw (11)<<height;
	dout<<setw (11)<<brchlen1;
	dout<<setw (7)<<num_child;
	dout<<setw (8)<<num_descndnt;
	dout<<setw(4)<<num_descndnt_interior;
	dout<<setw (6)<<rank<<"   ";
	//for (size_t i=0;i<descndnt.size();i++){
		//dout<<setw (1)<<descndnt[i];
	//}	
	dout<<setw(3)<<e_num;
	dout<<"    "<<clade;
	//dout<<setw(3)<<num_descndnt_interior;
}
