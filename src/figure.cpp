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

#include"figure.hpp"
//#include<exception>
#include <stdexcept>      // std::invalid_argument


Figure::Figure( int argc, char *argv[] ){
    this->init();
	for ( int argc_i = 1; argc_i < argc ; argc_i++ ){
		std::string argv_i( argv[argc_i] );
		if ( argv_i == "-label"  ){ this->check_option();   this->option = LABEL;  }
		if ( argv_i == "-branch" ){ this->check_option();   this->option = BRANCH; }
		if ( argv_i == "-dot"    ){ 
            this->check_method();   
            this->method = DOT;    
            this->figure_file_prefix = "dotfigure.dot";
            //this->figure_file_prefix = "dotfigure";
            this->figure_file_suffix = ".dot";
            }
		if ( argv_i == "-plot"   ){ 
            this->check_method();   
            this->method = LATEX;  
            this->figure_file_prefix = "texfigure.tex";
            this->figure_file_suffix = ".tex";
            }
        
		if (argv_i=="-dot_file" || argv_i=="-dotF"){
			this->check_method();   this->method = DOT; 
			this->figure_file_prefix = argv[argc_i+1];
            this->figure_file_suffix = ".dot";
			argc_i++;
		}
		

		if (argv_i=="-plot_file" || argv_i=="-plotF"){
			this->check_method();   this->method = LATEX; 
			this->figure_file_prefix = argv[argc_i+1];
			this->figure_file_suffix = ".tex";
            argc_i++;
		}
	}
    this->finalize();
}

void Figure::init(){
     this->option = PLOT_DEFAULT;
     this->method = NO_METHOD;
    }

void Figure::check_option(){
    if ( this->option != PLOT_DEFAULT ){
        throw std::invalid_argument ( " Plot option can either be \"-label\" or \"-branch\" " );
    }
}

void Figure::check_method(){
    if ( this->method != NO_METHOD ){
        throw std::invalid_argument ( " Method can either be LaTex (\"-plot\") or DOT (\"-dot\") " );
    }
}

/*! Check and remove files*/
void Figure::finalize(){
    size_t found =  this->figure_file_prefix.find( this->figure_file_suffix );
	if ( found != string::npos ){
        this->figure_file_prefix = this->figure_file_prefix.substr(0,found);
	}
    this->figure_file_name = this->figure_file_prefix + this->figure_file_suffix;
	
    ifstream tmp_file( this->figure_file_name.c_str() );
	if ( tmp_file.good() ) 	{  remove(figure_file_name.c_str()); 	}
}


/*! \brief Produce a tex file, which is used to draw the network 
 */
void Figure::plot_in_latex_file( ){
	ofstream latex_file;
    latex_file.open( this->figure_file_name.c_str(), ios::out | ios::app | ios::binary ); 
	latex_file << "\\documentclass[10pt]{article}\n";
	latex_file << "\\usepackage{tikz,graphics,graphicx,lscape,fullpage,multicol,setspace}\n \\singlespacing\n \\begin{document}\n ";	
	latex_file << "\\ifx\\du\\undefined\\newlength{\\du}\\fi\\setlength{\\du}{30\\unitlength}\n";
	latex_file << "\\begin{center}\n";
	latex_file.close();
    this->plot_in_latex_core();
    latex_file.open ( this->figure_file_name.c_str(), ios::out | ios::app | ios::binary ); 
	latex_file << "\\end{center}\n";
	latex_file << "\\end{document}\n";
	latex_file.close();
	
	string command = "pdflatex " + this->figure_file_name;
	int sys=system(command.c_str());
    std::clog << "Network figure generated in file: " + this->figure_file_name << endl; 
}

void Figure::plot_in_latex_core(){
	ofstream latex_file;
	latex_file.open (this->figure_file_name.c_str(), ios::out | ios::app | ios::binary); 	
	latex_file << "\\begin{tikzpicture}[thick]\n";
	valarray <int>  x_node = this->det_x_node ( );
	for (size_t node_i = 0; node_i < this->obj_net.Net_nodes.size();node_i++){
		string sp_node_label = this->obj_net.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		if (obj_net.Net_nodes[node_i].tip_bool){
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<obj_net.Net_nodes[node_i].rank<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		//latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<y_node[node_i]<<"\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
		else{
			latex_file<<"\\node at ("<<x_node[node_i]<<"\\du,"<<obj_net.Net_nodes[node_i].rank<<"\\du) [circle,fill=orange,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
	}

	for (size_t node_i=0; node_i < obj_net.Net_nodes.size()-1; node_i++ ){
		string sp_node_label=obj_net.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		string sp_node_parent1_label=obj_net.Net_nodes[node_i].parent1->label;
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!obj_net.Net_nodes[node_i].tip_bool){
			if (this->option == LABEL){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< obj_net.Net_nodes[node_i].e_num <<"};\n";
			}
			else{
				if (this->option == BRANCH){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< obj_net.Net_nodes[node_i].brchlen1 <<"};\n";
				}
				else{
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
				}	
			}
			if (obj_net.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=obj_net.Net_nodes[node_i].parent2->label;	
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (this->option == LABEL){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< obj_net.Net_nodes[node_i].e_num2 <<"};\n";
				}
				else{
					if (this->option == BRANCH){
						latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< obj_net.Net_nodes[node_i].brchlen2 <<"};\n";
					}
					else{
						latex_file<<"\\draw ("<<sp_node_label<<")-- (" << sp_node_parent2_label<<");\n";
					}		
				}
			}
		}
		else{
			if (this->option == BRANCH){
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<") node [midway,left]{"<< obj_net.Net_nodes[node_i].brchlen1 <<"};\n";
			}
			else{
				latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent1_label<<");\n";	
			}	
			if (obj_net.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=obj_net.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (this->option == BRANCH){
					latex_file<<"\\draw ("<<sp_node_label<<")-- ("<<sp_node_parent2_label<<") node [midway,left]{"<< obj_net.Net_nodes[node_i].brchlen2 <<"};\n";
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

/*! \brief Produce a dot file, which is used to draw the network, and compile the dot file to a pdf file.
 */
void Figure::plot_in_dot( ){
	ofstream dot_file;
	dot_file.open ( this->figure_file_name.c_str(), ios::out | ios::app | ios::binary ); 
	dot_file <<"graph G {\n rankdir=BT; ratio=compress;\n";//page="14,14"; determines the size of the ps output

	for ( size_t node_i = 0; node_i < this->obj_net.Net_nodes.size()-1; node_i++ ){    
		string sp_node_label = this->obj_net.Net_nodes[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		string sp_node_parent1_label=obj_net.Net_nodes[node_i].parent1->label;
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);
		if (!obj_net.Net_nodes[node_i].tip_bool){		
			if (this->option == LABEL){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< obj_net.Net_nodes[node_i].e_num <<"\"];\n";
			    }
			else if (this->option == BRANCH){
                dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< obj_net.Net_nodes[node_i].brchlen1 <<"\"];\n";	
				}
            else{
                dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
            }	
			
			if (obj_net.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=obj_net.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
				if (this->option == LABEL){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< obj_net.Net_nodes[node_i].e_num2 <<"\"];\n";
				    }
				else if (this->option == BRANCH){
                    dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< obj_net.Net_nodes[node_i].brchlen2 <<"\"];\n";	
					}
                else{
                    dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< obj_net.Net_nodes[node_i].e_num2 <<"\"];\n";
					}	
			}
		}
		else{
			if (this->option == BRANCH){
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<"[label=\""<< obj_net.Net_nodes[node_i].brchlen1 <<"\"];\n";	
			}
			else{
				dot_file<<sp_node_label<<" -- "<<sp_node_parent1_label<<";\n";//
			}	
			
			if (obj_net.Net_nodes[node_i].parent2){
				string sp_node_parent2_label=obj_net.Net_nodes[node_i].parent2->label;
				sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);					
				if (this->option == BRANCH){
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<"[label=\""<< obj_net.Net_nodes[node_i].brchlen2 <<"\"];\n";	
				}
				else{
					dot_file<<sp_node_label<<" -- "<<sp_node_parent2_label<<";\n";//<<"[label=\""<< obj_net.Net_nodes[node_i].e_num2 <<"\"];\n";
				}
			}
		}
	}
	
	if (obj_net.is_ultrametric){
		for (int rank_i=obj_net.Net_nodes.back().rank;rank_i>0;rank_i--){
			dot_file<<"{ rank=same; ";
			vector <int> x_node_dummy;
			vector <size_t> x_node_dummy_index;
			for (size_t node_i=0;node_i<obj_net.Net_nodes.size();node_i++){
				if (obj_net.Net_nodes[node_i].rank==rank_i){
					string sp_node_label=obj_net.Net_nodes[node_i].label;
					sp_node_label=rm_and_hash_sign(sp_node_label);
					dot_file<<sp_node_label<<" ";
				}	
			}
			dot_file<<"} ;\n";
		}
	}
	dot_file <<"}\n";
	dot_file.close();

	//string command = "dot -Tps2 " + this->figure_file_prefix + ".dot -o " + this->figure_file_prefix + ".eps";
	//int sys = system(command.c_str());
	//command = "dot -Tps2 " + this->figure_file_prefix + ".dot -o " + this->figure_file_prefix + ".ps";
	//sys = system(command.c_str());
    //command = "dot -Tpdf " + this->figure_file_prefix + ".dot -o " + this->figure_file_prefix + ".pdf";
	//sys=system(command.c_str());
    this->execute_dot ("ps2", ".eps");
    this->execute_dot ("ps2", ".ps");
    this->execute_dot ("pdf", ".pdf");
	std::clog << "Dot figure generated in file: " + this->figure_file_prefix + ".pdf" << endl;
}

void Figure::execute_dot(string method, string suffix){
	string command = "dot -T" + method + " " + this->figure_file_prefix + ".dot -o " + this->figure_file_prefix + suffix;
	(void)system( command.c_str() );
    }


/*! \brief When drawing network in .tex files, detemine the x coordinates of nodes
 */
valarray <int>  Figure::det_x_node ( ){
	valarray <int>  x_node ( this->obj_net.Net_nodes.size() );
	x_node[x_node.size()-1]=0;
	for ( int rank_i = this->obj_net.Net_nodes.back().rank; rank_i > 0; rank_i-- ){ 
		vector <int> x_node_dummy;
		vector <size_t> x_node_dummy_index;
		for ( size_t node_i = 0; node_i < this->obj_net.Net_nodes.size(); node_i++ ){
			if ( this->obj_net.Net_nodes[node_i].rank == rank_i){
				size_t n_child = this->obj_net.Net_nodes[node_i].child.size();
				int parent_x = x_node[node_i];
				int start_child_x = parent_x-floor(n_child/2);

                bool odd_num_child = (n_child % 2) == 1 ? true:false;
                
				if (odd_num_child){
					for (size_t child_i=0; child_i < n_child;child_i++){
						for (size_t node_j=0; node_j < this->obj_net.Net_nodes.size(); node_j++){
							//if ( this->obj_net.Net_nodes[node_j].label == this->obj_net.Net_nodes[node_i].child[child_i]->label){			
                            if ( &this->obj_net.Net_nodes[node_j] == this->obj_net.Net_nodes[node_i].child[child_i] ){
								if ( start_child_x == parent_x){
									x_node[node_j] = parent_x;
									start_child_x++;
								}
								else{
									x_node[node_j] = start_child_x;
								}
								start_child_x++;
							}
						}
					}
				}
				else{
					for (size_t child_i=0; child_i < n_child; child_i++){
						for ( size_t node_j = 0; node_j < this->obj_net.Net_nodes.size(); node_j++){
							//if (this->obj_net.Net_nodes[node_j].label == this->obj_net.Net_nodes[node_i].child[child_i]->label){
                            if ( &this->obj_net.Net_nodes[node_j] == this->obj_net.Net_nodes[node_i].child[child_i] ){
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
        //cout <<"stop"
		if (x_node_dummy.size() > 1){
			bool need_to_shift=true;
			while (need_to_shift){
				for (size_t x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (size_t x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
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
				for (size_t x_node_dummy_i=0;x_node_dummy_i<x_node_dummy.size();x_node_dummy_i++){
					int current_x_node_dummy=x_node_dummy[x_node_dummy_i];
					for (size_t x_node_dummy_j=x_node_dummy_i+1;x_node_dummy_j<x_node_dummy.size();x_node_dummy_j++){
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


void Figure::plot( string net_str ){
	this->obj_net = Net (net_str);
    //obj_net.print_all_node();
	if ( this->method == LATEX ){ this->plot_in_latex_file(); }
	if ( this->method == DOT   ){ this->plot_in_dot();		  }
	
	if ( this->option == BRANCH ){
		std::clog << std::endl << "Internal branches are labelled by post-order tree traversal." << std::endl;
	}
	if ( this->option == LABEL ){
		std::clog << std::endl << "Branch lengths are labelled." << std::endl;
	}
}
