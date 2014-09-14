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

bool Node::print_dout( bool is_Net ){
    dout << setw(12) << this;
	if ( is_Net ) dout << setw(6) << this->hybrid();
    if ( is_Net ) dout << setw(8) << descndnt_of_hybrid;
	dout << setw(5) << tip_bool;
    //if (this->parent1) dout << setw (11) << (parent1->label);
    //else dout << "           ";
    if (this->parent1) dout << setw (11) << (parent1);
    else dout << "           ";
	dout << setw (6) << this->height();
	dout << setw (12) << this->brchlen1();
    if (is_Net){
        if (this->parent2) dout << setw (11) << (parent2->label);
        else dout << "           ";
        dout<<setw (12) << this->brchlen2();
    }
	dout << setw (7) << this->child.size();
	dout << setw (8) << num_descndnt;
	dout << setw(4) << num_descndnt_interior;
	dout << setw(6) << this->rank() << "   ";
	//for (size_t i=0;i<descndnt.size();i++){
		//dout<<setw (1)<<descndnt[i];
	//}
	dout << setw(2)<<this->e_num();
	if ( is_Net ) dout << setw(3) << this->e_num2();
	dout << "    " << this->clade;
	//dout<<endl;
    return true;
}
