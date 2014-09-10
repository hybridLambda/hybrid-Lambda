
/*! \brief free the meomory */
void Net::clear(){
	tax_name.clear();
	NodeContainer.clear();
};

string Net::rewrite_internal_node_content( size_t i ){
    string new_node_content="(";
    for (size_t child_i = 0; child_i < this->NodeContainer[i].child.size(); child_i++ ){
        if ( this->NodeContainer[i].child[child_i]->node_content == this->NodeContainer[i].child[child_i]->label ) {
            new_node_content += this->NodeContainer[i].child[child_i]->label+":" + to_string ( this->NodeContainer[i].child[child_i]->brchlen1() ) ;
        }
        else {
            string brchlen_str2;
            bool new_hybrid_node=false;
            for ( size_t node_ii=0; node_ii < i; node_ii++){
                for ( size_t node_ii_child_i = 0; node_ii_child_i < this->NodeContainer[node_ii].child.size(); node_ii_child_i++ ){
                    if ( this->NodeContainer[node_ii].child[node_ii_child_i]->node_content == this->NodeContainer[i].child[child_i]->node_content){
                        new_hybrid_node=true;
                        brchlen_str2 = to_string(this->NodeContainer[i].child[child_i]->brchlen2() );
                    break;}
                }
                if (new_hybrid_node){break;}
            }
            if ( !new_hybrid_node ) new_node_content += this->NodeContainer[i].child[child_i]->node_content;
            //new_node_content += this->NodeContainer[i].child[child_i]->label+":" + brchlen_str2;
            new_node_content += new_hybrid_node ? this->NodeContainer[i].child[child_i]->label+":" + brchlen_str2 : 
                                                  this->NodeContainer[i].child[child_i]->node_content + this->NodeContainer[i].child[child_i]->label+":" + brchlen_str2 ;
        }
        if ( child_i < this->NodeContainer[i].child.size() - 1 ) new_node_content += ",";
    }
    new_node_content += ")";
    return new_node_content;
}
