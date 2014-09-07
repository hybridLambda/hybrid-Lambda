void print_all_child(Node *parent /*! pointer to the parent node*/){
    cout << parent->label << " has " << parent->num_child << " kids" << endl;
    for (int i_num_child=0;i_num_child<=parent->num_child-1;i_num_child++){
        cout<<parent->child[i_num_child]->label<<endl;
    }
}


void print_parent(Node *child /*! pointer to the child node*/){
	cout<<child->label<<" has parents "<<endl;
	cout<<child->parent1->label<<" and "<<child->parent2->label <<endl;
}
