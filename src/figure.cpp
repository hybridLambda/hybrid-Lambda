/* 
 * 
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012, 2013 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-Lambda 
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

#include"figure.hpp"

figure::param::param(){
	 this->plot_bool=false;
	 dot_bool=false;
	 plot_option=0;
	 plot_label=false;
	 plot_branch=false;
	 this->tex_fig_name="texfigure.tex";
	 this->dot_fig_name="dotfigure.dot";
}

figure::param::param(int argc, char *argv[]){
	 plot_bool=false;
	 dot_bool=false;
	 
	 plot_option=0;
	 plot_label=false;
	 plot_branch=false;
	 tex_fig_name="texfigure.tex";
	 dot_fig_name="dotfigure.dot";
	for (int argc_i=1; argc_i < argc ;argc_i++){
		
		std::string argv_i(argv[argc_i]);
		if (argv_i=="-label"){
			plot_label=true;
		}

		if (argv_i=="-branch"){
			plot_branch=true;
		}		

		if (argv_i=="-dot"){
			dot_bool=true;
		}
		if (argv_i=="-dot_file" || argv_i=="-dotF"){
			dot_bool=true;
			dot_fig_name=argv[argc_i+1];
			argc_i++;
		}
		check_and_remove(dot_fig_name.c_str());
		
		if (argv_i=="-plot"){
			plot_bool=true;
		}
		if (argv_i=="-plot_file" || argv_i=="-plotF"){
			plot_bool=true;
			tex_fig_name=argv[argc_i+1];
			argc_i++;

		}
		check_and_remove(tex_fig_name.c_str());
	}
	set_plot_option_();
}

void figure::param::plot(string net_str){
	Net net_dummy(net_str);
	if (plot_bool){
		plot_in_latex_file_(net_dummy);	
	}
	if (dot_bool){
		plot_in_dot_(net_dummy);			
	}
	
	if (plot_option==1){
		std::cout<<std::endl<<"Internal branches are labelled by post-order tree traversal."<<std::endl;
	}
	if (plot_option==2){
		std::cout<<std::endl<<"Branch lengths are labelled."<<std::endl;
	}
}

/*! \brief Set plot options */
/*! \return int '0' -- do not label the branches; '1' -- enumerate the branches by  postorder tree traversal; '2' -- label the branch lengths of the network*/ 
//int figure::param::set_plot_option_(bool plot_label,bool plot_branch){
void figure::param::set_plot_option_(){
	if (plot_label){
		plot_option=1;
		
	}
	else{
		if (plot_branch){
			plot_option=2;		
		}
	}
}

/*! \brief When drawing network in .tex files, detemine the x coordinates of nodes
 */
valarray <int>  figure::param::det_x_node (Net net_dummy){
	valarray <int>  x_node (net_dummy.Net_nodes.size());
	x_node[x_node.size()-1]=0;
	
	for (int rank_i=net_dummy.Net_nodes.back().rank;rank_i>0;rank_i--){
		vector <int> x_node_dummy;
		vector <unsigned int> x_node_dummy_index;
		for (unsigned int node_i=0;node_i<net_dummy.Net_nodes.size();node_i++){
			if (net_dummy.Net_nodes[node_i].rank==rank_i){
				unsigned int n_child=net_dummy.Net_nodes[node_i].child.size();
				int parent_x=x_node[node_i];
				int start_child_x=parent_x-floor(n_child/2);
				bool odd_num_child=false;
				if ((n_child % 2) == 1){
					odd_num_child=true;
				}
				if (odd_num_child){
					for (unsigned int child_i=0;child_i<n_child;child_i++){
						for (unsigned int node_j=0;node_j<net_dummy.Net_nodes.size();node_j++){
							if (net_dummy.Net_nodes[node_j].label==net_dummy.Net_nodes[node_i].child[child_i]->label){			
								if (start_child_x==parent_x){
									x_node[node_j]=parent_x;
									start_child_x++;
								}
								else{
									x_node[node_j]=start_child_x;
								}
								start_child_x++;
							}
						}
					}
				}
				else{
					for (unsigned int child_i=0;child_i<n_child;child_i++){
						for (unsigned int node_j=0;node_j<net_dummy.Net_nodes.size();node_j++){
							if (net_dummy.Net_nodes[node_j].label==net_dummy.Net_nodes[node_i].child[child_i]->label){
								if (start_child_x==parent_x){										
									start_child_x++;
								}
								//	x_node[node_j]=start_child_x*(parent_y-child_y)+parent_x;
								x_node[node_j]=start_child_x;
								//y_node[node_j]=parent_y-1;
								start_child_x++;
							}
						}
					}
					
				}
				x_node_dummy.push_back(x_node[node_i]);
				x_node_dummy_index.push_back(node_i);
			}
		}
		if (x_node_dummy.size()>1){
			bool need_to_shift=true;
			while (need_to_shift){
				for (unsigned int x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (unsigned int x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
						if (current_x_node_dummy==x_node_dummy[x_node_dummy_j]){
							if (x_node_dummy[x_node_dummy_j]>0){
								x_node_dummy[x_node_dummy_j]++;
								x_node[x_node_dummy_index[x_node_dummy_j]]++;
							}
							else{
								x_node_dummy[x_node_dummy_j]--;
								x_node[x_node_dummy_index[x_node_dummy_j]]--;
							}
						}
					}
				}
				need_to_shift=false;
				for (unsigned int x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (unsigned int x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
						if (current_x_node_dummy==x_node_dummy[x_node_dummy_j]){
							need_to_shift=true;
							break;
						}		
					}
					if (need_to_shift){
						break;
					}
				}
			}
		}
	}
	return x_node;
	
}




/*! \brief Produce a dot file, which is used to draw the network, and compile the dot file to a pdf file.
 */
void figure::param::plot_in_dot_(//const char* file_name /*! Name for the figure file */,
	Net net_dummy
	){
	string file_name_no_dot(dot_fig_name.c_str());
	string file_name_with_dot(dot_fig_name.c_str());
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
	sys=system(command.c_str());
	string appending_log_str="Dot figure generated in file: "+file_name_no_dot+".pdf";
}



/*! \brief Core function of drawing a network in .tex files. 
 */
 void figure::param::plot_in_latex_(const char* file_name /*! Name for the figure file */ , 
	Net net_dummy
	){
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
 void figure::param::plot_in_latex_file_(//const char* file_name /*! Name for the figure file */ , 
	Net net_dummy
	){
	ofstream latex_file;
	string file_name_no_dot(tex_fig_name.c_str());
	string file_name_with_dot(tex_fig_name.c_str());
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
	//plot_in_latex(file_name_with_dot.c_str(), net_dummy,plot_option);	
	plot_in_latex_(file_name_with_dot.c_str(),net_dummy);	
	//		plot_in_latex(file_name, net_str,plot_option);	
	latex_file.open (file_name_with_dot.c_str(), ios::out | ios::app | ios::binary); 
	latex_file <<"\\end{center}\n";
	latex_file <<"\\end{document}\n";
	latex_file.close();
	
	string command="pdflatex "+file_name_no_dot+".tex";
	int sys=system(command.c_str());
	string appending_log_str="Network figure generated in file: "+file_name_no_dot+".pdf";

}

