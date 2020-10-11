/*
 * hybrid-Lambda is used to simulate gene trees given species network under
 * coalescent process.
 *
 * Copyright (C) 2010 -- 2015 Sha (Joe) Zhu
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

#include"node.hpp"

Node::Node(){
    label="";
    node_content="";
    //num_child=0;
    num_descndnt=0;
    num_descndnt_interior=0;
    parent1=NULL;
    parent2=NULL;
    this->brchlen1_ = 0.0;
    this->brchlen2_ = 0.0;
    this->rank_     = 0;
    this->e_num_    = 0;
    this->e_num2_   = 0;
    //visit=0;
    descndnt_of_hybrid=false;
    tip_bool=false;
    //clade=" ";
    this->height_ = 1.0/0.0;
    //prob_to_hybrid_left=1.0;
    this->visited_ = false;
    this->set_label1_starts_at(0);
    this->set_label2_starts_at(0);
}


void Node::print( bool is_Net ){
    cout << setw(12) << label;
    if ( is_Net ) cout << setw(6) << this->hybrid();
    if ( is_Net ) cout << setw(8) << descndnt_of_hybrid;
    cout << setw(5) << tip_bool;
    if (this->parent1) cout << setw (11) << (parent1->label);
    else cout << "           ";
    cout << setw (12) << this->height();
    cout << setw (12) << this->brchlen1();
    if (is_Net){
        if (this->parent2) cout << setw (11) << (parent2->label);
        else cout << "           ";
        cout<<setw (12) << this->brchlen2();
    }
    cout << setw (7) << this->child.size();
    cout << setw (8) << num_descndnt;
    cout << setw(4) << num_descndnt_interior;
    cout << setw(6) << this->rank() << "   ";
    //for (size_t i=0;i<descndnt.size();i++){
        //cout<<setw (1)<<descndnt[i];
    //}
    cout << setw(2)<<this->e_num();
    if ( is_Net ) cout << setw(3) << this->e_num2();
    cout << "    " << this->clade;
    //cout<<endl;
}

/*! \brief Add child node to parent node */
void Node::add_child( Node *child_node, /*! pointer to the child node*/
                      size_t adding_to_parent){
    this->child.push_back(child_node);

    if ( (adding_to_parent == 2) | (child_node->parent1 != NULL)){
        child_node->parent2 = (this);
        //child_node->hybrid = true;
    } else {
        assert(adding_to_parent == 1);
        child_node->parent1 = (this);
    }
}


/*! \brief Rank network node from the bottom.
 *
 * Child node has lower rank than the parent node. Tip nodes have rank one, the root node has the highest rank
 */
void Node::CalculateRank(){
    if ( this->tip_bool ) this->rank_ = 1;
    else {
        size_t child_max_rank = 0;
        for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++ ){
            this->child[ith_child]->CalculateRank();
            child_max_rank = max( child_max_rank, this->child[ith_child]->rank() );
        }
        this->rank_ = child_max_rank + 1;
    }
}


bool Node::find_descndnt ( string &name, NAMETYPE type ){
    if ( this->tip_bool ) {
        string tmp = ( type == TAXA ) ? this->name : this->label;
        return ( tmp == name ) ? true : false;
    }
    else {
        bool descndnt_found = false;
        for (size_t i = 0; i < this->child.size(); i++ ){
            descndnt_found = this->child[i]->find_descndnt( name, type );
            if ( descndnt_found ) break;
        }
        return descndnt_found;
    }
}


/*! \brief label a node if its a tip node */
void Node::find_tip(){
    if ( this->child.size() == 0) this->tip_bool = true;
    else {
        for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++ ){
            (*this->child[ith_child]).find_tip();
            //this->child[ith_child]->find_tip();
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////// consider removed

/*! \brief Label a node if its a descendant of a hybrid node */
void Node::find_hybrid_descndnt(){
    if ( this->tip_bool ) return;
    for ( size_t ith_child = 0; ith_child < this->child.size(); ith_child++){
        if ( this->hybrid() || this->descndnt_of_hybrid ) this->child[ith_child]->descndnt_of_hybrid = true;
        this->child[ith_child]->find_hybrid_descndnt();
    }
}
