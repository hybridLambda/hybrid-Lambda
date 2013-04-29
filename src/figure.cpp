//figure.cpp


/*! \brief Produce a dot file, which is used to draw the network, and compile the dot file to a pdf file.
 */
void plot_in_dot(const char* file_name /*! Name for the figure file */,
	Net net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/){
	string file_name_no_dot(file_name);
	string file_name_with_dot(file_name);
	string ddot(".dot");
	
	size_t found=file_name_no_dot.find(ddot);
	if (found!=string::npos){
		file_name_no_dot=file_name_no_dot.substr(0,found);
	}
	else{
		file_name_with_dot=file_name_with_dot+ddot;
	}

	//Net net_dummy(net_str);
	ofstream dot_file;
	check_and_remove(file_name_with_dot.c_str());
	dot_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
			
	dot_file <<"graph G {\n rankdir=BT; ratio=compress;\n";//page="14,14"; determines the size of the ps output
	//valarray <int>  x_node=det_x_node (net_dummy);

	for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size()-1;node_i++){
		string sp_node_label=net_dummy.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		string sp_node_parent1_label=net_dummy.Net_nodes[node_i].parent1->label;
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!net_dummy.Net_nodes[node_i].tip_bool){		
			if (plot_option==1){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num <<"\"];\n";
			}
			else{
				if (plot_option==2){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen1 <<"\"];\n";	
				}
				else{
					dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
				}	
			}
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==1){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num2 <<"\"];\n";
				}
				else{
					if (plot_option==2){
						dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen2 <<"\"];\n";	
					}
					else{
						dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num2 <<"\"];\n";
					}	
				}	
				
			}
		}
		else{
			if (plot_option==2){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen1 <<"\"];\n";	
			}
			else{
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
			}	
			
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);					
				if (plot_option==2){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< net_dummy.Net_nodes[node_i].brchlen2 <<"\"];\n";	
				}
				else{
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< net_dummy.Net_nodes[node_i].e_num2 <<"\"];\n";
				}
			}
		}
	}
	
	if (net_dummy.is_ultrametric){
		for (int rank_i=net_dummy.Net_nodes.back().rank;rank_i>0;rank_i--){
			dot_file<<"{ rank=same; ";
			vector <int> x_node_dummy;
			vector <unsigned int> x_node_dummy_index;
			for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size();node_i++){
				if (net_dummy.Net_nodes[node_i].rank==rank_i){
					string sp_node_label=net_dummy.Net_nodes[node_i].label;
					sp_node_label=rm_and_hash_sign(sp_node_label);
					dot_file<<sp_node_label<<" ";
				}	
			}
			dot_file<<"} ;\n";
		}
	}
	dot_file <<"}\n";
	dot_file.close();

	string command="dot -Tps "+file_name_no_dot+".dot -o "+file_name_no_dot+".ps";
	int sys=system(command.c_str());
	command="convert "+file_name_no_dot+".ps -resize 100\% "+file_name_no_dot+".pdf";
	//cout<<command<<endl;
	sys=system(command.c_str());

	string appending_log_str="Dot figure generated in file: "+file_name_no_dot+".pdf";
	appending_log_file(appending_log_str);
}



/*! \brief Core function of drawing a network in .tex files. 
 */
void plot_in_latex(const char* file_name /*! Name for the figure file */ , 
	Net net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/
	){
	//Net net_dummy(net_str);
	ofstream latex_file;
	latex_file.open (file_name, ios::out | ios::app | ios::binary); 	
	latex_file <<"\\begin{tikzpicture}[thick]\n";
	valarray <int>  x_node=det_x_node (net_dummy);
	for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size();node_i++){
		string sp_node_label=net_dummy.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		if (net_dummy.Net_nodes[node_i].tip_bool){
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<net_dummy.Net_nodes[node_i].rank<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		//latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<y_node[node_i]<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
		else{
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<net_dummy.Net_nodes[node_i].rank<<"\\du) [circle,fill=orange,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
	}

	for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size()-1;node_i++){
		string sp_node_label=net_dummy.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		string sp_node_parent1_label=net_dummy.Net_nodes[node_i].parent1->label;
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!net_dummy.Net_nodes[node_i].tip_bool){
			if (plot_option==1){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].e_num <<"};\n";
			}
			else{
				if (plot_option==2){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen1 <<"};\n";
				}
				else{
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
				}	
			}
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;	
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==1){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].e_num2 <<"};\n";
				}
				else{
					if (plot_option==2){
						latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen2 <<"};\n";
					}
					else{
						latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<");\n";
					}		
				}
			}
		}
		else{
			if (plot_option==2){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen1 <<"};\n";
			}
			else{
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
			}	
			if (net_dummy.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=net_dummy.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (plot_option==2){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< net_dummy.Net_nodes[node_i].brchlen2 <<"};\n";
				}
				else{
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<");\n";
				}		
			}
		}
	}
	latex_file <<"\\end{tikzpicture}\n\n";
	latex_file.close();	
}



/*! \brief Produce a tex file, which is used to draw the network 
 */
void plot_in_latex_file(const char* file_name /*! Name for the figure file */ , 
	Net net_dummy,
// string net_str /*! Input network written in extended newick form */,
	int plot_option /*! '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/ 
	){
	ofstream latex_file;
	string file_name_no_dot(file_name);
	string file_name_with_dot(file_name);
	string dtex(".tex");
	
	size_t found=file_name_no_dot.find(dtex);
	if (found!=string::npos){
		file_name_no_dot=file_name_no_dot.substr(0,found);
	}
	else{
		file_name_with_dot=file_name_with_dot+dtex;
	}
	
	latex_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
	latex_file <<"\\documentclass[10pt]{article}\n";
	latex_file <<"\\usepackage{tikz,graphics,graphicx,lscape,fullpage,multicol,setspace}\n \\singlespacing\n \\begin{document}\n ";	
	latex_file<<"\\ifx\\du\\undefined\\newlength{\\du}\\fi\\setlength{\\du}{30\\unitlength}\n";
	latex_file <<"\\begin{center}\n";
	latex_file.close();
	plot_in_latex(file_name_with_dot.c_str(), net_dummy,plot_option);	
	//		plot_in_latex(file_name, net_str,plot_option);	
	latex_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
	latex_file <<"\\end{center}\n";
	latex_file <<"\\end{document}\n";
	latex_file.close();
	
	string command="pdflatex "+file_name_no_dot+".tex";
	int sys=system(command.c_str());

	string appending_log_str="Network figure generated in file: "+file_name_no_dot+".pdf";
	appending_log_file(appending_log_str);

}

