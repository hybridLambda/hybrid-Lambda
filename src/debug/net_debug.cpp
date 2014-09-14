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

#include "../net.hpp"

bool Tree::print_all_node_dout(){

    if ( is_Net ) dout<<"           label  hybrid hyb_des non-tp parent1  height brchln1 parent2 brchln2 #child #dsndnt #id rank   e_num   Clade "<<endl;
    else dout<<"       label non-tp   parent  height brchln #child #dsndnt #id rank e_num   Clade "<<endl;
    for (size_t i=0;i<NodeContainer.size();i++){
        //for (size_t j=0;j<descndnt[i].size();j++) dout<<setw(3)<<descndnt[i][j];

        assert( this->NodeContainer[i].print_dout( this->is_Net_() ) );

        dout<<"  ";
        
        for (size_t j=0;j<samples_below[i].size();j++) dout<<samples_below[i][j];        
        dout<<endl;
    }
    return true;
}
