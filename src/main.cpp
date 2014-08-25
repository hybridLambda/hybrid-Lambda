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

/*! \file main.cpp
 *  \brief Main function of hybrid-Lambda */

#include "hybridLambda.hpp"
#include "fst.hpp"

int main(int argc, char *argv[]){

	if ( argc == 1 ) print_help(); 	//else, proceed

    try {
	    HybridLambda run_hybridLambda ( argc, argv );
        double Fst;

        run_hybridLambda.HybridLambda_core( );
        run_hybridLambda.create_site_data_dir(); // segregating site data were generated				

        run_hybridLambda.extract_tmrca ();
        run_hybridLambda.extract_bl ();
        run_hybridLambda.extract_firstcoal();
        run_hybridLambda.extract_frequency();


/// need to work         
//if (run_hybridLambda.fst_bool){
////sim::param sim_para(argc, argv);	
//Net coal_unit_net(sim_para.net_str);
//double tau = coal_unit_net.Net_nodes[0].brchlen1;

////Net para_net(sim_para.para_string);
////double lambdaA  = para_net.Net_nodes[0].brchlen1;
////double lambdaAB = para_net.Net_nodes.back().brchlen1;
////cout<<"tau = " << tau <<endl;
////cout << "( 1 - exp( -tau ) )= "<<( 1 - exp( -tau ) )<<endl;
//Fst = ( 1 - exp( -tau ) ) * ( tau / ( 1 + tau ) );
//cout << "Expected[Fst] = " << Fst << endl;
//}

		//}			
    }
    catch (const exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}

void print_example(){
	cout<<"Examples:"<<endl;
	cout<<""<<endl;	
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 3 -seed 2 -o example1"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp_nt1_para -o example2 -num 2 -mu 0.00003 -sim_mut_unit -sim_num_mut"<<endl;
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 100 -pop 25000 -sim_num_gener"<<endl;
	cout<<"hybrid-Lambda -spng '(A:50000,B:50000)r;' -pop '(A:50000,B:50000)r:40000;'"<<endl;
	cout<<"hybrid-Lambda -spcu '((((A:1.1,B:1.1):2.1,a:2.2):1.1,13D:.2):.3,4:.3);' -S 2 4 3 6 5"<<endl;
	cout<<"hybrid-Lambda -spcu '(A:1,B:1)r;' -mm '(A:1.9,B:.2)r:2;' -S 3 4"<<endl;
	cout<<"hybrid-Lambda -spcu trees/7_tax_sp_nt1_para -dot -branch"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp1.tre -num 1000 -o GENE -f"<<endl;	
	//cout<<"hybrid-Lambda -spcu trees/4_tax_sp1.tre -num 1000 -o GENE_TREE_FILE -fF FRENQUENCY_FILE"<<endl;	
	//cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 1000 -o GENE -fF OUTPUT"<<endl;	
	cout<<"hybrid-Lambda -gt GENE_coal_unit -f "<<endl;	
    cout<<"hybrid-Lambda -mt GENE_num_mut -seg "<<endl;	
	cout<<"hybrid-Lambda -spcu '(A:5,B:5)r;' -mono -num 100 -mm .1 -S 4 4"<<endl;
	cout<<endl;
}

void print_option(){
	cout<<setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<endl;
	cout<<setw(20)<<"-spcu STR"<<"  --  "<<"Input the species network/tree string through command line or a file."<<endl;
	cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in coalescent unit."<<endl;
	cout<<setw(20)<<"-spng STR"<<"  --  "<<"Input the species network/tree string through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in number of generation."<<endl;
    cout<<setw(20)<<"-gt STR"<<"  --  "<<"Specify the FILE of trees to analyse tree topology frequencies."<<endl;
    cout<<setw(20)<<"-mt STR"<<"  --  "<<"Specify the FILE of trees to analyse tree topology frequencies."<<endl;
	cout<<setw(20)<<"-pop STR/FLOAT"<<"  --  "<<"Population sizes are defined by a single numerical constant, "<<endl;
	cout<<setw(26)<<" "<<"or a string which specifies the population size on each branch. "<<endl;
	cout<<setw(26)<<" "<<"The string can be input through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"By default, population size 10,000 is used."<<endl;
	cout<<setw(20)<<"-mm STR/FLOAT"<<"  --  "<<"Multiple merger parameters are defined by a single numerical constant, "<<endl;
	cout<<setw(26)<<" "<<"or a string which speifies the parameter on each branch. "<<endl;
	cout<<setw(26)<<" "<<"The string can be input through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"By default, Kingman coalescent is used."<<endl;
	cout<<setw(20)<<"-S INT INT ..."<<"  --  "<<"Specify the number of samples for each taxon."<<endl;
	cout<<setw(20)<<"-num INT"<<"  --  "<<"The number of gene trees will be simulated."<<endl;
	cout<<setw(20)<<"-seed INT"<<"  --  "<<"User define random SEED"<<endl;
	cout<<setw(20)<<"-mu FLOAT"<<"  --  "<<"User defined constant mutation rate MU. By default mutation rate 0.00005 is used."<<endl;
	cout<<setw(20)<<"-o STR [option]"<<"  --  "<<"Specify the file name prefix for simulated gene trees. \"GENE_TREE\" by default"<<endl;
	//cout<<"     By default, gene tree branch lengths are in coalescent unit "<<endl;
	cout<<setw(20)<<"-sim_mut_unit"<<"  --  "<<"Convert the simulated gene tree branch lengths to mutation unit."<<endl;
	cout<<setw(20)<<"-sim_num_gener"<<"  --  "<<"Convert the simulated gene tree branch lengths to number of generations."<<endl;
	cout<<setw(20)<<"-sim_num_mut"<<"  --  "<<"Simulate numbers of mutations on each branch of simulated gene trees."<<endl;
	cout<<setw(20)<<"-sim_Si_num"<<"  --  "<<"Generate the file out table, which includes the number of segregating"<<endl;
	cout<<setw(26)<<" "<<"sites and the total branch length of the gene tree in coalescent unit."<<endl;
	cout<<setw(20)<<"-f"<<"  --  "<<"Generate a topology frequency table of a set of input trees or simulated gene trees."<<endl;
	cout<<setw(26)<<" "<<"Frequency table is saved in file freq out by default."<<endl;
	//cout<<setw(20)<<"-fF FILE"<<"  --  "<<"The topology frequency table will be saved in the FILE."<<endl;	
	cout<<setw(20)<<"-mono"<<"  --  "<<"Generate a frequency table of monophyletic, paraphyletic and polyphyletic trees. "<<endl;
	cout<<setw(20)<<"-plot/-dot [option]"<<"  --  "<<"Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree)."<<endl;
	cout<<setw(20)<<"      -branch"<<"  --  "<<"Branch lengths will be labelled in the figure."<<endl;
	//cout<<setw(20)<<"-plotF/-dotF FILE"<<"  --  "<<"Generated figure will be saved in FILE."<<endl;			
	cout<<endl;	
}

/*! \brief hybrid-Lambda help file*/
void print_help(){
	cout<<endl;
	cout<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<"*                      hybrid-Lambda beta 0.4                   *"<<endl;
	cout<<"*                         Author: Joe ZHU                       *"<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<endl<<endl;
	print_option();
	print_example();
    exit (EXIT_SUCCESS);
}
