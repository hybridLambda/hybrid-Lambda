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

#include"net.hpp"

/*! \brief Construct Net object from a (extended) Newick string */
Tree::Tree(string old_string /*! input (extended) newick form string */){
    if (old_string.size()==0){
        descndnt.clear();
        tax_name.clear();
        NodeContainer.clear();
        return;
    }

    this->init();
    this->check_Parenthesis(old_string);
    this->check_labeled( old_string );
    // check & sign, this should be illigal for hybrid-Lambda,

    vector<string> labels;
    vector<string> node_contents;
    vector<string> brchlens;
    vector<size_t> label_starts_at;
    vector<size_t> content_starts_at;
    this->initialize_nodes_label_and_content(labels, node_contents, brchlens, label_starts_at, content_starts_at);
    this->initialize_NodeContainer(labels, node_contents, brchlens, label_starts_at, content_starts_at);
    this->extract_tax_and_tip_names();

    this->connect_graph();
    this->NodeContainer.back().find_tip();
    this->NodeContainer.back().find_hybrid_descndnt();
    this->NodeContainer.back().CalculateRank();
    this->max_rank = NodeContainer.back().rank();
    this->enumerate_internal_branch( this->NodeContainer.back() );
    this->init_descendant();
    this->init_node_clade();
    //this->rewrite_descendant();
    this->check_isNet();
    this->check_isUltrametric();
    //dout<<"Net constructed"<<endl;
}


void Tree::initialize_nodes_label_and_content(vector<string> & labels,
                                              vector<string> & node_contents,
                                              vector<string> & brchlens,
                                              vector<size_t> & label_starts_at,
                                              vector<size_t> & content_starts_at) {
    size_t found_bl = net_str.find(':');
    for (size_t i_str_len=1;i_str_len<net_str.size();){
        if (net_str[i_str_len]=='e' && (net_str[i_str_len+1]=='-' || net_str[i_str_len+1]=='+')){
            i_str_len++;
        }
        else{
            if ( start_of_tax_name(net_str,i_str_len) ){
                size_t str_start_index = i_str_len;
                string label = extract_label(net_str,i_str_len);
                labels.push_back(label);
                label_starts_at.push_back(str_start_index);

                string node_content;
                if ( net_str[str_start_index-1]==')' ){
                    size_t rev_dummy_i = Parenthesis_balance_index_backwards( net_str, str_start_index-1 );
                    size_t substr_len = str_start_index-rev_dummy_i;
                    node_content = net_str.substr(rev_dummy_i, substr_len );
                    content_starts_at.push_back(rev_dummy_i);
                }
                else {
                    node_content=label;
                    content_starts_at.push_back(str_start_index);
                }
                i_str_len += label.size();

                node_contents.push_back(node_content);
                string brchlen;
                if ( found_bl != string::npos ){
                    size_t found=min(min(net_str.find(",",i_str_len+1),net_str.find(")",i_str_len+1)),net_str.size());
                    brchlen = net_str.substr(i_str_len+1,found-i_str_len-1);
                }
                found_bl = net_str.find(":", found_bl+1);
                brchlens.push_back(brchlen);
            }
            else {
                i_str_len++;
            }
        }
    }
}


void Tree::initialize_NodeContainer(vector<string> & labels,
                                    vector<string> & node_contents,
                                    vector<string> & brchlens,
                                    vector<size_t> & label_starts_at,
                                    vector<size_t> & content_starts_at){
    //int label_counter = brchlens.size();
    for ( size_t new_i_label=0 ; new_i_label < brchlens.size(); new_i_label++ ){
        Node empty_node;
        NodeContainer.push_back(empty_node);
        NodeContainer[new_i_label].label = labels[new_i_label];
        NodeContainer[new_i_label].node_content = node_contents[new_i_label];
        NodeContainer[new_i_label].set_brchlen1( strtod(brchlens[new_i_label].c_str(), NULL) );
        //cout << new_i_label << "  " << labels[new_i_label] << "  " << brchlens[new_i_label] << endl;
        NodeContainer[new_i_label].set_label1_starts_at(label_starts_at[new_i_label]);
        NodeContainer[new_i_label].set_node_content_starts_at(content_starts_at[new_i_label]);
    }

    for ( size_t i = 1; i < NodeContainer.size()-1; i++ ){
        size_t j;
        for ( j = i+1; j < NodeContainer.size()-1; j++ ){
            if ( NodeContainer[j].label==NodeContainer[i].label ){
                if ( NodeContainer[j].node_content[0] == '(' ){
                    NodeContainer[i].node_content = NodeContainer[j].node_content;
                }
                NodeContainer[i].set_brchlen2 ( NodeContainer[j].brchlen1() );
                NodeContainer[i].set_label2_starts_at ( NodeContainer[j].label1_starts_at() );
                NodeContainer[i].set_node_content_starts_at ( NodeContainer[j].node_content_starts_at() );
                break;
            }
        }
        if ( NodeContainer[j].label == NodeContainer[i].label ) NodeContainer.erase(NodeContainer.begin()+j);
    }

}


void Tree::init_descendant(){
    for ( size_t i = 0; i < NodeContainer.size(); i++){
        valarray <int> descndnt_dummy(0,tax_name.size());
        descndnt.push_back(descndnt_dummy);
        valarray <int> samples_below_dummy(0,tip_name.size());
        samples_below.push_back(samples_below_dummy);
        for ( size_t tax_name_i = 0; tax_name_i < tax_name.size(); tax_name_i++ ) descndnt[i][tax_name_i] = this->NodeContainer[i].find_descndnt( tax_name[tax_name_i], TAXA) ? 1:0;
        for ( size_t tip_name_i = 0; tip_name_i < tip_name.size(); tip_name_i++ ) samples_below[i][tip_name_i] = this->NodeContainer[i].find_descndnt( tip_name[tip_name_i], TIP) ? 1:0;
        this->NodeContainer[i].num_descndnt = descndnt[i].sum();
    }

    for ( size_t i = 0; i < NodeContainer.size(); i++){
        for (size_t j = 0; j < NodeContainer.size(); j++){
            if ( i == j ) continue;

            valarray <int> descndnt_diff=(descndnt[i]-descndnt[j]);
            if (descndnt_diff.min() >= 0 && NodeContainer[i].rank() > NodeContainer[j].rank() && NodeContainer[j].rank() >= 2){
                this->NodeContainer[i].num_descndnt_interior += 1 ;
                this->NodeContainer[i].descndnt_interior_node.push_back( &this->NodeContainer[j] );
            }
        }
    }
}


void Tree::init_node_clade(){
    for ( size_t i = 0; i < NodeContainer.size(); i++ ){
        //if ( this->descndnt[i].sum() == 0 ) break;

        this->NodeContainer[i].clade.clear();
        for ( size_t tax_name_i = 0; tax_name_i < tax_name.size(); tax_name_i++ ){
            if ( descndnt[i][tax_name_i] != 1) continue;

            NodeContainer[i].clade = ( NodeContainer[i].clade.size() == 0 ) ? tax_name[tax_name_i]:
                                                                            NodeContainer[i].clade + tax_name[tax_name_i];
            NodeContainer[i].clade.push_back('&');
        }
        NodeContainer[i].clade.erase(NodeContainer[i].clade.size()-1,1);
    }
}

string Tree::extract_label(string in_str, size_t i){
    size_t j=end_of_label_or_bl(in_str, i);
    //cout<<"i="<<i<<", j="<<j<<endl;
    return in_str.substr(i,j+1-i);
}



void Tree::extract_tax_and_tip_names(){
    for (size_t i=0;i<NodeContainer.size();i++){
        if(NodeContainer[i].label != NodeContainer[i].node_content) continue;
        if ( NodeContainer[i].label.find("_") > 0 ){
            //multi_label_bool=true;
            NodeContainer[i].name = NodeContainer[i].label.substr(0,NodeContainer[i].label.find("_"));
            //cout<<NodeContainer[i].name<<endl;
            bool new_tax_bool=true;
            for ( size_t tax_i = 0; tax_i < tax_name.size(); tax_i++ ){
                if (tax_name[tax_i]==NodeContainer[i].name){
                    new_tax_bool=false;
                    break;
                }
            }
            if ( new_tax_bool ){
                tax_name.push_back(NodeContainer[i].name);
            }
            //cout<<tax_name.back()<<endl;
        }
        else{
            tax_name.push_back(NodeContainer[i].label);
        }
        tip_name.push_back(NodeContainer[i].label);
    }
    sort(tax_name.begin(), tax_name.end());
    sort(tip_name.begin(), tip_name.end());
    //cout << " tax_name.size() = "<<tax_name.size()<<endl;
    //cout << " tip_name.size() = "<<tax_name.size()<<endl;
}


void Tree::connect_graph(){
    for ( size_t i = 0; i < NodeContainer.size(); i++ ){
        if ( NodeContainer[i].node_content[0] != '(' ) continue;

        // cout <<  NodeContainer[i].label << " " <<  NodeContainer[i].node_content <<endl;
        char child_node1[NodeContainer[i].node_content.length()];
        for ( size_t i_content_len = 1; i_content_len < NodeContainer[i].node_content.length(); ){
            if (NodeContainer[i].node_content[i_content_len]=='(' ||  start_of_tax_name(NodeContainer[i].node_content,i_content_len) ){
                size_t j_content_len = (NodeContainer[i].node_content[i_content_len] == '(') ? Parenthesis_balance_index_forwards( NodeContainer[i].node_content, i_content_len ) + 1:
                                                                                               i_content_len;
                int child1_node_content_i = 0;
                for ( ; j_content_len < NodeContainer[i].node_content.length(); j_content_len++){
                    child_node1[child1_node_content_i] = NodeContainer[i].node_content[j_content_len];
                    char stop = NodeContainer[i].node_content[j_content_len+1];
                    if ( stop == ',' || stop == ')' || stop == ':'){
                        child_node1[child1_node_content_i+1]='\0';
                        break;
                    }
                    child1_node_content_i++;
                }

                string child_node1_str = child_node1;
                i_content_len = j_content_len + 2;
                for ( size_t j = 0; j < NodeContainer.size(); j++){

                    size_t adding_to_parent = 1;

                    if (child_node1_str == NodeContainer[j].label) {
                        // cout<< NodeContainer[i].label << " " << NodeContainer[i].node_content_starts_at() + child1_node_content_i << " " << NodeContainer[j].label << "  " <<  NodeContainer[j].label1_starts_at() <<endl;
                        // Checking node content label index ...
                        if (NodeContainer[j].label2_starts_at() == 0){
                            assert(adding_to_parent = 1);
                        } else {
                            if (NodeContainer[i].node_content_starts_at() > NodeContainer[j].label1_starts_at()){
                                adding_to_parent = 2;
                            }
                        }


                        NodeContainer[i].add_child( &NodeContainer[j], adding_to_parent );
                    }
                }
            }
            else { i_content_len++;}
        }
    }
}


void Tree::check_labeled( string in_str ){
    bool labeled_bool=true;
    for ( size_t i = 0; i < in_str.size(); i++ ){
        if ( in_str[i] == ')' && i == end_of_label_or_bl(in_str, i) ){
            labeled_bool = false;
            break;
        }
    }
    this->net_str = labeled_bool ? in_str:label_interior_node(in_str);
}


void Tree::check_isUltrametric(){
    vector <int> remaining_node( NodeContainer.size(), 0 );
    for ( size_t node_i = 0; node_i < NodeContainer.size(); node_i++ ){
        remaining_node[node_i] = node_i;
    }
    size_t rank_i = 1;
    size_t remaining_node_i=0;
    while ( remaining_node.size() > 0 ){
        int node_i = remaining_node[remaining_node_i];
        if ( NodeContainer[node_i].rank() == rank_i ){
            if (rank_i == 1) NodeContainer[node_i].path_time.push_back(0.0);
            else{
                for (size_t child_i = 0; child_i < NodeContainer[node_i].child.size(); child_i++ ){

                    double current_child_time = (NodeContainer[node_i].child[child_i]->parent1->label==NodeContainer[node_i].label)?
                                                NodeContainer[node_i].child[child_i]->brchlen1():
                                                NodeContainer[node_i].child[child_i]->brchlen2();
                    for (size_t child_i_time_i=0;child_i_time_i<NodeContainer[node_i].child[child_i]->path_time.size();child_i_time_i++){
                        NodeContainer[node_i].path_time.push_back(current_child_time+NodeContainer[node_i].child[child_i]->path_time[child_i_time_i]);
                    }
                }
            }
            remaining_node.erase(remaining_node.begin()+remaining_node_i);
        }
        else{
            remaining_node_i++;
        }

        if ( remaining_node_i == remaining_node.size()-1 ){
            rank_i++;
            remaining_node_i=0;
        }
    }

    for (size_t node_i=0;node_i<NodeContainer.size();node_i++){
        for (size_t path_time_i=0;path_time_i<NodeContainer[node_i].path_time.size();path_time_i++){
            if (pow((NodeContainer[node_i].path_time[path_time_i]-NodeContainer[node_i].path_time[0]),2)>0.000001){
                this->is_ultrametric = false;
                break;
            }
        }
        NodeContainer[node_i].set_height( NodeContainer[node_i].path_time[0] );
    }
}


void Tree::check_isNet(){ //false stands for tree, true stands for net_work
    for (size_t i = 0; i < this->NodeContainer.size(); i++){
        if ( !this->NodeContainer[i].parent2 ) continue;
        this->is_Net = true;
        return;
    }
}


void Tree::print_all_node(){
    if ( this->is_Net ) cout<<"           label  hybrid hyb_des non-tp parent1  abs_t brchln1 parent2 brchln2 #child #dsndnt #id rank   e_num   Clade "<<endl;
    else cout<<"            label non-tp   parent        abs_t brchln #child #dsndnt #id rank e_num   Clade "<<endl;
    for (size_t i = 0; i < this->NodeContainer.size(); i++ ){
        for (size_t j = 0; j < this->descndnt[i].size(); j++ ) {cout<<setw(3)<<this->descndnt[i][j];}
        this->NodeContainer[i].print( this->is_Net_() );
        cout<<"  ";
        for (size_t j=0;j<this->samples_below[i].size();j++) {cout<<this->samples_below[i][j]; }
        cout<<endl;
    }
}


/*! \brief Label interior node if the interior nodes of the tree string are not labeled */
string Tree::label_interior_node(string in_str /*!< input newick form string */){
    vector <string> in_str_partition;
    int interior_node_counter = 0;
    int sub_str_start_index = 0;
    size_t i = in_str.find(')');
    while ( i<in_str.size() ){
        interior_node_counter++;
        string current_string;
        size_t found_next_bracket = min(in_str.find(")", sub_str_start_index),in_str.size());
        current_string = in_str.substr(sub_str_start_index, found_next_bracket - sub_str_start_index +1);
        if ( in_str[i+1] == ';' || i == (in_str.size()-1) ){
            current_string += "root";
        }
        else {
            current_string += "Int_" + to_string(static_cast<int>(interior_node_counter));
            sub_str_start_index = i+1;
        }
        in_str_partition.push_back(current_string);

        i = in_str.find( ")", i+1 );
    }
    string out_str;
    for ( size_t i = 0; i < in_str_partition.size(); i++ )     out_str += in_str_partition[i];
    return out_str;
}


size_t Tree::first_coal_rank(){
    size_t min_rank = NodeContainer.back().rank();
    for (size_t i = 0 ; i < NodeContainer.size(); i++){
        if ( NodeContainer[i].tip_bool ) continue;
        min_rank = ( NodeContainer[i].rank() < min_rank ) ?  NodeContainer[i].rank() : min_rank ;
    }
    return min_rank;
}


size_t Tree::first_coal_index (){
    size_t min_rank = this->first_coal_rank();
    size_t dummy_index = this->NodeContainer.size()-1;
    double min_coal_time = this->NodeContainer[dummy_index].height();
    for (size_t i = 0 ; i < NodeContainer.size(); i++){
        if ( this->NodeContainer[i].rank() == min_rank &&  this->NodeContainer[i].height() < min_coal_time ){
            dummy_index = i;
            min_coal_time = this->NodeContainer[dummy_index].height();
        }
    }
    return dummy_index;
}


/*! \brief enumerate the internal branches */
void Tree::enumerate_internal_branch( Node & node ) {
    if ( node.tip_bool ) return;

    if ( node.visited() ){
        this->current_enum_ ++;
        node.set_enum2( current_enum_ );
        }
    else{
        for ( size_t i = 0; i < node.child.size(); i++ ){
            this->enumerate_internal_branch( *node.child[i] );
        }
        node.set_visited( true );
        this->current_enum_ ++;
        node.set_enum( current_enum_ );
    }
}


/*! \brief Identify if its the start of the taxon name in a newick string, should be replaced by using (isalpha() || isdigit())  */
bool Tree::start_of_tax_name( string in_str, size_t i ){
    //bool start_bool = false;
    //if ( (in_str[i]!='(' && in_str[i-1]=='(') || (in_str[i-1]==',' && in_str[i]!='(') || ( (in_str[i-1]==')') && ( in_str[i]!=')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) ) {
        //start_bool=true;
    //}
    //return     start_bool;
    if      (  in_str[i-1] == '('  &&   in_str[i] != '(' ) return true;
    else if (  in_str[i-1] == ','  &&   in_str[i] != '(' ) return true;
    else if ( (in_str[i-1] == ')') && ( in_str[i] != ')' || in_str[i]!=':' || in_str[i]!=',' || in_str[i]!=';' ) ) return true;
    else return false;
}


size_t Tree::Parenthesis_balance_index_backwards( string &in_str, size_t i ){
    size_t j = i;
    int num_b = 0;
    for ( ; j > 0 ; j-- ){
        if      ( in_str[j] == '(' ) num_b--;
        else if ( in_str[j] == ')' ) num_b++;
        else continue;
        if ( num_b == 0 ) break;
    }
    return j;
}


size_t Tree::Parenthesis_balance_index_forwards( string &in_str, size_t i ){
    size_t j = i;
    int num_b = 0;
    for ( ; j < in_str.size(); j++ ){
        if      ( in_str[j] == '(' ) num_b++;
        else if ( in_str[j] == ')' ) num_b--;
        else continue;
        if ( num_b == 0 ) break;
    }
    return j;
}


/*! \brief Checking Parenthesis of a (extended) Newick string */
void Tree::check_Parenthesis( string &in_str ){
    int num_b = 0;
    for ( size_t i = 0; i < in_str.size(); i++){
        if      (in_str[i] == '(') num_b++;
        else if (in_str[i] == ')') num_b--;
        else continue;
    }
    if ( num_b != 0 ) throw std::invalid_argument(in_str + "Parenthesis not balanced!" );
}



/*! \brief rewrite node content of nodes */
void Tree::rewrite_node_content(){
    int highest_i = 0;
    for ( size_t i = 0; i < this->NodeContainer.size(); i++ ){
        if ( this->NodeContainer[i].num_descndnt > this->NodeContainer[highest_i].num_descndnt ){ highest_i = i;}
    }

    this->NodeContainer[highest_i].CalculateRank();

    for ( size_t rank_i = 1; rank_i <= this->NodeContainer.back().rank(); rank_i++){
        for ( size_t i = 0 ; i < this->NodeContainer.size(); i++ ){
            if ( this->NodeContainer[i].rank() != rank_i ) continue;

            this->NodeContainer[i].node_content = ( this->NodeContainer[i].rank() == 1 ) ?
                                                    this->NodeContainer[i].label :
                                                    this->rewrite_internal_node_content( i );
        }
    }
}


string Tree::rewrite_internal_node_content( size_t i ){
    string new_node_content="(";
    for (size_t child_i = 0; child_i < this->NodeContainer[i].child.size(); child_i++ ){
        string brchlen_str1 = to_string (static_cast<double>(this->NodeContainer[i].child[child_i]->brchlen1() ));
        if ( this->NodeContainer[i].child[child_i]->node_content == this->NodeContainer[i].child[child_i]->label ) {
            //new_node_content += this->NodeContainer[i].child[child_i]->label+":" + to_string ( this->NodeContainer[i].child[child_i]->brchlen1() ) ;
            new_node_content += this->NodeContainer[i].child[child_i]->label+":" +  brchlen_str1;
        }
        else {
            bool new_hybrid_node=false;
            string brchlen_str2;
            for ( size_t node_ii=0; node_ii < i; node_ii++){
                for ( size_t node_ii_child_i = 0; node_ii_child_i < this->NodeContainer[node_ii].child.size(); node_ii_child_i++ ){
                    if ( this->NodeContainer[node_ii].child[node_ii_child_i]->node_content == this->NodeContainer[i].child[child_i]->node_content){
                        new_hybrid_node=true;
                        brchlen_str2 = to_string(static_cast<double>(this->NodeContainer[i].child[child_i]->brchlen2() ));
                    break;}
                }
                if (new_hybrid_node){break;}
            }
            new_node_content += new_hybrid_node ? this->NodeContainer[i].child[child_i]->label+":" + brchlen_str2 :
                                                  this->NodeContainer[i].child[child_i]->node_content + this->NodeContainer[i].child[child_i]->label+":" +  brchlen_str1;
        }
        if ( child_i < this->NodeContainer[i].child.size() - 1 ) new_node_content += ",";
    }
    new_node_content += ")";
    return new_node_content;
}

/////////////////////////////////////////// consider for removal

void Tree::rewrite_descendant(){    //check for coaleased tips(& sign in the tips)
    bool rewrite_descndnt=false;
    for (size_t i=0;i<NodeContainer.size();i++){
        if (NodeContainer[i].tip_bool ){
            for (size_t i_str=0;i_str<NodeContainer[i].clade.size();i_str++){
                if (NodeContainer[i].clade[i_str]=='&'){
                    rewrite_descndnt=true;
                    break;
                }
            }
        }
        if (rewrite_descndnt){
            break;
        }
    }

    if ( !rewrite_descndnt ) return;

    tax_name.clear();
    int tax_name_start=0;
    int tax_name_length=0;
    for (size_t new_i_str=0;new_i_str<NodeContainer.back().clade.size();new_i_str++){
        tax_name_length++;
        if (NodeContainer.back().clade[new_i_str]=='&'){
            tax_name_length--;
            tax_name.push_back(NodeContainer.back().clade.substr(tax_name_start,tax_name_length));
            tax_name_start=new_i_str+1;
            tax_name_length=0;
        }
        if (new_i_str==NodeContainer.back().clade.size()-1){
            tax_name.push_back(NodeContainer.back().clade.substr(tax_name_start,tax_name_length));
        }
    }
    sort(tax_name.begin(), tax_name.end());
    descndnt.clear();

    for (size_t i=0;i<NodeContainer.size();i++){
        vector <string> contained_tips;
        valarray <int> re_initial_descndnt(0,tax_name.size());
        int tax_name_start=0;
        int tax_name_length=0;
        for (size_t new_i_str=0;new_i_str<NodeContainer[i].clade.size();new_i_str++){
            tax_name_length++;
            if (NodeContainer.back().clade[new_i_str]=='&'){
                tax_name_length--;
                contained_tips.push_back(NodeContainer[i].clade.substr(tax_name_start,tax_name_length));
                tax_name_start=new_i_str+1;
                tax_name_length=0;
            }
            if (new_i_str==NodeContainer[i].clade.size()-1){
                contained_tips.push_back(NodeContainer[i].clade.substr(tax_name_start,tax_name_length));
            }
        }
        for (size_t tax_i=0;tax_i<tax_name.size();tax_i++){
            for (size_t contained_tax_i=0;contained_tax_i<contained_tips.size();contained_tax_i++){
                if (tax_name[tax_i]==contained_tips[contained_tax_i]){
                    //descndnt[i][tax_i]=1;
                    re_initial_descndnt[tax_i]=1;
                }
            }
        }
        descndnt.push_back(re_initial_descndnt);
    }
    //this->rewrite_node_clade();
    this->init_node_clade();
}

string Tree::print_newick( Node * node ){
    string tree_str;
    if ( node->tip_bool ){
        tree_str = node->label ;
    }
    else {
        tree_str = "(";
        for ( size_t i = 0 ; i < node->child.size() ; i++ ){
            tree_str += print_newick ( node->child[i] ) + ":" + to_string (static_cast<double>(node->child[i]->brchlen1() ));
            if ( i < node->child.size()-1 ) tree_str += ",";
        }
        tree_str += ")";
    }
    return tree_str;
}


/*! \brief Remove interior nodes label of a string */
string remove_interior_label(string in_str/*!< input newick form string */){
    string out_str;
    out_str=in_str;

    size_t found_bracket=out_str.find(')');
    while ( found_bracket<out_str.size() ){
        if ( isalpha(out_str[found_bracket+1]) || isdigit(out_str[found_bracket+1]) ){
            size_t char_j = end_of_label_or_bl( out_str, found_bracket+1 );
            out_str.erase(out_str.begin()+found_bracket+1, out_str.begin()+char_j+1);
        }
        found_bracket = out_str.find( ")",found_bracket+1 );
    }
    return out_str;
}


size_t end_of_label_or_bl( string &in_str, size_t i ){
    //size_t j ;
    //for ( j = i; j < in_str.size(); j++){
        //char stop=in_str[j+1];
        //if ( stop == ',' || stop==')' || stop==':' || stop==';'){
            //break;
        //}
    //}
    //return j;
    size_t j = i;
    for ( ; j < in_str.size(); j++){
        if      ( in_str[j+1] == ',' )    break;
        else if ( in_str[j+1] == ')' )    break;
        else if ( in_str[j+1] == ':' )    break;
        else if ( in_str[j+1] == ';' )    break;
        else continue;
    }
    return j;
}


void readNextStringto( string &readto , int& argc_i, int argc_, char * const* argv_ ){
    argc_i++;
    if (argc_i >= argc_) throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv_[argc_i-1]);
    readto = std::string(argv_[argc_i]);
    if ( readto[0] == '-' ) throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv_[argc_i-1]);
}
