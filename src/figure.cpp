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
#include <stdexcept>      // std::invalid_argument

Figure::Figure( int argc, char* const *argv ):
    argc_(argc), argv_(argv){
    this->init();
	while( argc_i < argc_ ){
		std::string argv_i( argv_[argc_i] );
		if ( argv_i == "-label"  ){ this->check_option();   this->option = LABEL;  }
		if ( argv_i == "-branch" ){ this->check_option();   this->option = BRANCH; }        
		if ( argv_i == "-dot"    ){ 
            this->check_method();  
            this->initialize_method ( DOT, ".dot" );
            this->figure_file_prefix = "dotfigure";
        }
		if ( argv_i == "-plot"   ){ 
            this->check_method();   
            this->initialize_method ( LATEX, ".tex" );
            this->figure_file_prefix = "texfigure";
        }
		if (argv_i == "-dot_file" || argv_i == "-dotF"){
			this->check_method();
            this->initialize_method ( DOT, ".dot" );
            //this->read_prefix();
            readNextStringto( this->figure_file_prefix , this->argc_i, this->argc_,  this->argv_ );
		}
		if (argv_i == "-plot_file" || argv_i == "-plotF"){
			this->check_method();
            this->initialize_method ( LATEX, ".tex" );
            //this->read_prefix();
            readNextStringto( this->figure_file_prefix , this->argc_i, this->argc_,  this->argv_ );
		}
        argc_i++;
	}
    this->finalize();
}

void Figure::initialize_method( FIGURE_PROGRAM program, string suffix){
    this->method = program;  
    this->figure_file_suffix = suffix;
    }


void Figure::check_option(){
    if ( this->option != PLOT_DEFAULT ){
        cerr << "Too many figure options!"<<endl;
        throw std::invalid_argument ( " Plot option can either be \"-label\" or \"-branch\" " );
    }
}

void Figure::check_method(){
    if ( this->method != NO_METHOD ){
        cerr << "Which method (figure) do you mean?"<<endl;
        throw std::invalid_argument ( " Method can either be LaTex (\"-plot\" or \"-plotF\") or DOT (\"-dot\" or \"-dotF\") " );
    }
}

void Figure::init(){
     this->option = PLOT_DEFAULT;
     this->method = NO_METHOD;
     this->argc_i = 1;
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
void Figure::plot_in_latex( ){
    figure_ofstream.open( this->figure_file_name.c_str(), ios::out | ios::app | ios::binary ); 
	figure_ofstream << "\\documentclass[10pt]{article}\n";
	figure_ofstream << "\\usepackage{tikz,graphics,graphicx,lscape,fullpage,multicol,setspace}\n \\singlespacing\n \\begin{document}\n ";	
	figure_ofstream << "\\ifx\\du\\undefined\\newlength{\\du}\\fi\\setlength{\\du}{30\\unitlength}\n";
	figure_ofstream << "\\begin{center}\n";
    figure_ofstream << "\\begin{tikzpicture}[thick]\n";
    this->det_x_node ( );
	for (size_t node_i = 0; node_i < this->obj_net.NodeContainer.size();node_i++){
		string sp_node_label = this->obj_net.NodeContainer[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		if (obj_net.NodeContainer[node_i].tip_bool){
			figure_ofstream<<"\\node at ("<<x_node[node_i]<<"\\du,"<<obj_net.NodeContainer[node_i].rank() << "\\du) [circle,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
		else{
			figure_ofstream<<"\\node at ("<<x_node[node_i]<<"\\du,"<<obj_net.NodeContainer[node_i].rank() << "\\du) [circle,fill=orange,draw] ("<<sp_node_label <<") {$"<<sp_node_label <<"$};\n";
		}
	}
    this->plot_core();	    
   	figure_ofstream <<"\\end{tikzpicture}\n\n";
	figure_ofstream << "\\end{center}\n";
	figure_ofstream << "\\end{document}\n";
	figure_ofstream.close();
	
	string command = "pdflatex " + this->figure_file_name;
	int sys=system(command.c_str());
    std::clog << "Network figure generated in file: " + this->figure_file_name << endl; 
}

void Figure::edge_entry(string from, string to, size_t label, double bl, bool tip){
    if ( this->option == LABEL && tip ){ 
        figure_ofstream << ( ( this->method == DOT )  ? "" : "\\draw (" )
                        << from 
                        << ( ( this->method == DOT )  ? " -- " : ")--(" )
                        << to 
                        << ( ( this->method == DOT )  ? "[label=\"" : ") node [midway,left]{" )
                        << label
                        << ( ( this->method == DOT )  ? "\"];\n" : "};\n" );
    }
    else if (this->option == BRANCH){
        figure_ofstream << ( ( this->method == DOT )  ? "" : "\\draw (" )
                        << from 
                        << ( ( this->method == DOT )  ? " -- " : ")--(" )
                        << to 
                        << ( ( this->method == DOT )  ? "[label=\"" : ") node [midway,left]{" )
                        << bl
                        << ( ( this->method == DOT )  ? "\"];\n" : "};\n" );
    }
    else{
        figure_ofstream << ( ( this->method == DOT ) ? "" : "\\draw (" )
                        << from
                        << ( ( this->method == DOT ) ? " -- " : ")--(" )
                        << to 
                        << ( ( this->method == DOT ) ? ";\n" : ");\n" );
    }	
}


void Figure::plot_core(){
   	for ( size_t node_i = 0; node_i < this->obj_net.NodeContainer.size()-1; node_i++ ){    
		string sp_node_label = this->obj_net.NodeContainer[node_i].label;
		sp_node_label=rm_and_hash_sign(sp_node_label);
		string sp_node_parent1_label=obj_net.NodeContainer[node_i].parent1->label;
		sp_node_parent1_label=rm_and_hash_sign(sp_node_parent1_label);

        this->edge_entry(sp_node_label, sp_node_parent1_label, obj_net.NodeContainer[node_i].e_num(), obj_net.NodeContainer[node_i].brchlen1(), !obj_net.NodeContainer[node_i].tip_bool);
        if (obj_net.NodeContainer[node_i].parent2){
            string sp_node_parent2_label=obj_net.NodeContainer[node_i].parent2->label;
            sp_node_parent2_label=rm_and_hash_sign(sp_node_parent2_label);
            this->edge_entry( sp_node_label, sp_node_parent2_label, obj_net.NodeContainer[node_i].e_num2(), obj_net.NodeContainer[node_i].brchlen2(), !obj_net.NodeContainer[node_i].tip_bool);
        }
	}
}

/*! \brief Produce a dot file, which is used to draw the network, and compile the dot file to a pdf file.
 */
void Figure::plot_in_dot( ){
    figure_ofstream.open ( this->figure_file_name.c_str(), ios::out | ios::app | ios::binary ); 
	figure_ofstream <<"graph G {\n rankdir=BT; ratio=compress;\n";//page="14,14"; determines the size of the ps output

    this->plot_core();	

	if (this->obj_net.is_ultrametric){
		for (size_t rank_i = this->obj_net.NodeContainer.back().rank(); rank_i > 0; rank_i--){
			figure_ofstream<<"{ rank=same; ";
			for (size_t node_i=0;node_i<obj_net.NodeContainer.size();node_i++){
				if (obj_net.NodeContainer[node_i].rank() == rank_i){
					string sp_node_label = obj_net.NodeContainer[node_i].label;
					sp_node_label = rm_and_hash_sign(sp_node_label);
					figure_ofstream << sp_node_label << " ";
				}	
			}
			figure_ofstream<<"} ;\n";
		}
	}
	figure_ofstream << "}\n";
	figure_ofstream.close();

    this->execute_dot ("ps2", ".eps");
    this->execute_dot ("ps2", ".ps");
    this->execute_dot ("pdf", ".pdf");
	std::clog << "Dot figure generated in file: " + this->figure_file_prefix + ".pdf" << endl;
}

void Figure::execute_dot(string method, string suffix){
	string command = "dot -T" + method + " " + this->figure_file_prefix + ".dot -o " + this->figure_file_prefix + suffix;
	int sys = system( command.c_str() );
    }


/*! \brief When drawing network in .tex files, detemine the x coordinates of nodes
 */
void  Figure::det_x_node ( ){
    this->x_node = valarray <int> (this->obj_net.NodeContainer.size());
	x_node[x_node.size()-1] = 0; //root x-axis value is zero

	for ( size_t rank_i = this->obj_net.NodeContainer.back().rank()-1; rank_i > 0; rank_i-- ){ // start from the children of the root
		this->x_node_tmp.clear();
		this->x_node_tmp_index.clear();
        
		for ( size_t node_i = 0; node_i < this->obj_net.NodeContainer.size(); node_i++ ){
			if ( this->obj_net.NodeContainer[node_i].rank() == rank_i){
				size_t n_child = this->obj_net.NodeContainer[node_i].child.size();
				int parent_x = x_node[node_i];
				int start_child_x = parent_x - floor( n_child / 2 );

                bool odd_num_child = (n_child % 2) == 1 ? true:false;
                dout << "parent x = " << parent_x <<" " << "start_child_x = "<<start_child_x <<" ";
                for ( size_t child_i = 0; child_i < n_child; child_i++ ){                    
                    dout << " child_"<<child_i << " x = " ;
                    for ( size_t node_j = 0; node_j < this->obj_net.NodeContainer.size(); node_j++ ){                    
                        if ( node_j == this->obj_net.NodeContainer[node_i].child[child_i]->node_index ){
                            if ( !odd_num_child && start_child_x == parent_x ){
                                start_child_x++;
                            }
                            x_node[node_j] = start_child_x;
                            start_child_x++;
                            dout << x_node[node_j] ;
                        }
                    }
                    
                }

				this->x_node_tmp.push_back(x_node[node_i]);
				this->x_node_tmp_index.push_back(node_i);
			}
		}
        this->x_node_shift();
        dout << endl;
	}
}

void Figure::x_node_shift(){
    if ( x_node_tmp.size() < 2 ){
    //cout <<"  no shift"<<endl;
        return;
    }
    //cout <<"  Need shift"<<endl;
    bool need_to_shift=true;
    while ( need_to_shift ){
        for (size_t x_node_tmp_i=0;x_node_tmp_i<x_node_tmp.size();x_node_tmp_i++){
            int current_x_node_tmp=x_node_tmp[x_node_tmp_i];
            for (size_t x_node_tmp_j=x_node_tmp_i+1;x_node_tmp_j<x_node_tmp.size();x_node_tmp_j++){
                if (current_x_node_tmp == x_node_tmp[x_node_tmp_j]){
                    if ( x_node_tmp[x_node_tmp_j] > 0 ){
                        x_node_tmp[x_node_tmp_j]++;
                        x_node[x_node_tmp_index[x_node_tmp_j]]++;
                    }
                    else{
                        x_node_tmp[x_node_tmp_j]--;
                        x_node[x_node_tmp_index[x_node_tmp_j]]--;
                    }
                }
            }
        }
        need_to_shift=false;
        for (size_t x_node_tmp_i=0;x_node_tmp_i<x_node_tmp.size();x_node_tmp_i++){
            int current_x_node_tmp=x_node_tmp[x_node_tmp_i];
            for (size_t x_node_tmp_j=x_node_tmp_i+1;x_node_tmp_j<x_node_tmp.size();x_node_tmp_j++){
                if (current_x_node_tmp==x_node_tmp[x_node_tmp_j]){
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


void Figure::plot( string net_str ){
    if ( net_str.size() == 0 ) { throw std::invalid_argument ("Population structure is undefined!!!");}
	this->obj_net = Net (net_str);
	if ( this->method == LATEX ){ this->plot_in_latex(); }
	if ( this->method == DOT   ){ this->plot_in_dot();   }
	
	if ( this->option == BRANCH ){
		std::clog << "Internal branches are labelled by post-order tree traversal." << std::endl;
	}
	if ( this->option == LABEL ){
		std::clog << "Branch lengths are labelled." << std::endl;
	}
}
