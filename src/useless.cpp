//useless.cpp


//AC::AC(){
	//pAC=0.0;
	//}
class CAC{
	public:
	vector < unsigned int > alive_lineages; 
	vector < class AC* > child_AC;
	
		unsigned int branch_index;
	double pCAC;
		CAC (){
		//pAC=0.0;
		pCAC=0.0;
	}
};
 
class AC{
	public:
	vector < unsigned int > alive_lineages; 
	vector < unsigned int > tip_starting_alive_lineages;
// 	vector < unsigned int > starting_alive_lineages; should be replaced by
	//vector < AC *> starting_AC;

	unsigned int branch_index;
	//vector < int > num_of_lineages_coal_child_AC;
	//vector < AC* > child_AC;
	vector < class CAC * > child_CAC;
	
	//vector < vector < AC* > > child_AC; /*! \todo !!! change to this, instead of vector < AC* > child_AC, some of AC may have multiple ways to be formed */
	double pAC;

	AC (){
		pAC=0.0;
		//pCAC=0.0;
	}
};






/*! \brief Compute factorial of a \return double a! */
double factorial (double a){
	if (a > 1){
		return (a * factorial (a-1));}
	else{
		return (1);}
}


/*! \brief Compute a permutations of n \return double */
double n_permu_a (double n, double a){
	if (a>1){
		return (n*n_permu_a(n-1,a-1));
	}
	else{
		if (a==1){
			return (n);
		}
		else{
			return (1);
		}
	}
}

/*! \brief Compute n choose k \return double */
double n_choose_k(double n, double k){
	if (k<(n/2)){
		return (n_choose_k(n,n-k));}
	else{
		return (n_permu_a(n,k)/factorial(k));}
}

/*! \brief Compute factorial of a \return int a! */
int factorial_int (int a){
	if (a > 1){
		return (a * factorial_int (a-1));}
	else{
		return (1);}
}

/*! \brief Compute a permutations of n \return int */
int n_permu_a_int (int n, int a){
	if (a>1){
		return (n*n_permu_a_int(n-1,a-1));}
	else{
		if (a==1){
			return (n);}
		else{
			return (1);}
	}
}

/*! \brief Compute n choose k \return int */
int n_choose_k_int(int n, int k){
	if (k<(n/2)){
		return (n_choose_k_int(n,n-k));}
	else{
		return (n_permu_a_int(n,k)/factorial_int(k));}
}


///*! \brief hybrid-Lambda help file*/
//void print_help(){
	//cout<<"*****************************************************************"<<endl;
	//cout<<"*			hybrid-Lambda beta 0.1			*"<<endl;
	//cout<<"*			  Author: Joe ZHU			*"<<endl;
	//cout<<"*****************************************************************"<<endl;
	//cout<<endl<<endl;
	//cout<<"-h or -help -- Help menu"<<endl;
	//cout<<"-num NUMBER -- Define the number of the gene tree will be simulated"<<endl;
	//cout<<"By default, hybrid-Lambda simulates gene trees under Kingmman Coalescent process, with using "<<endl;
	//cout<<"-mm -- it allows to have multi merger process"<<endl;
	//cout<<"-mm ALPHA -- 1 < ALPHA < 2 "<<endl;
	//cout<<"-mm PSI -- 0 < PSI < 1 "<<endl;
	//cout<<"-mm SPECIES_NETWORK_FILE_para -- Define the location of the SPECIES_NETWORK_FILE_para branch length are coalescent parameters in the branch"<<endl;
	//cout<<"-spcu SPECIES_NETWORK_FILE_coal_unit -- Define the location of the SPECIES_NETWORK_FILE branch length in coalescent unit"<<endl;
	//cout<<"-spng SPECIES_NETWORK_FILE_num_gener -- Define the location of the SPECIES_NETWORK_FILE branch length indicate number of generations in the branch"<<endl;
	//cout<<"-pop SPECIES_NETWORK_FILE_pop_size -- Define the location of the SPECIES_NETWORK_FILE branch length indicates population size in the branch"<<endl;
	//cout<<"-pop -- User define constant population size"<<endl;
	//cout<<"-mu -- User define mutation rate"<<endl;
	//cout<<"-gt GENE_TREE_FILE -- Define the location of the GENE_TREE_FILE"<<endl;
	//cout<<"-gF GENE_TREE_FILE -- Define the out name of file that gene trees will be saved, \"GENE_TREE\" by default"<<endl;
	//cout<<"     By default, gene tree branch lengths are in coalescent unit "<<endl;
	//cout<<"     -sim_mut_unit, gene tree branch lengths are in mutation unit "<<endl;
	//cout<<"     -sim_num_gener, gene tree branch lengths are in number of generations "<<endl;
	//cout<<"     -sim_num_mut, gene tree branch lengths are in number of mutations "<<endl;
	//cout<<"     -sim_Si_num, number of segregating sites "<<endl;
	//cout<<"-f -- Count the frequency of the simulated trees, output is saved in \"freq_out\" by default"<<endl;
	//cout<<"-fF FRENQUENCY_FILE -- Define the name of the file that count the frequency of the simulated trees"<<endl;
	//cout<<"-S NUMBER_OF_LINEAGE_ENTERING_1 NUMBER_OF_LINEAGE_ENTERING_1 ... -- Specify number of lineage entering each taxon "<<endl;
	//cout<<"-mono -- Give frequency of topology of monophyly, paraphyly and polyphyly"<<endl;
	//cout<<"-seed -- User define random seed"<<endl;
	//cout<<"-seg -- To produce segregating site data"<<endl;
	
	
	////cout<<"-debug -- To generate a debug file \"debug_file\""<<endl;
			
	//cout<<endl;
	//cout<<"Examples:"<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 "<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -dot "<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -dotF figure.dot "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plot "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plotF figure.tex "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plot -branch"<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -plotF -label figure.tex "<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000"<<endl;	
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -num 100 -S 100 100"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mm 1.1 -num 100 -S 100 100"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mm 0.1 -num 100 -S 100 100"<<endl;
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -gF GENE_TREE_FILE"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -f"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -fF FRENQUENCY_FILE"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -gF GENE_TREE_FILE -f"<<endl;	
	//cout<<"hybrid-Lambda -spcu 4_tax_sp1 -num 1000 -gF GENE_TREE_FILE -fF FRENQUENCY_FILE"<<endl;	
	//cout<<"hybrid-Lambda -gt GENE_TREE_FILE -f "<<endl;	
	//cout<<"hybrid-Lambda -gt GENE_TREE_FILE -fF FRENQUENCY_FILE"<<endl;	
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mono -num 100 -mm .1 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mono -num 100 -seed 2 -mm .1 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -para 2_tax_sp1_para2 -mono -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mm 2_tax_sp1_para2 -seg -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mm 2_tax_sp1_para2 -seg -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_size 10000 -para 2_tax_sp1_para2 -seg -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mu 0.0002 -para 2_tax_sp1_para2 -seg -num 100 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_string 2_tax_sp1_pop1 -mu 0.00002  -seg -num 1000 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spcu 2_tax_sp1 -mu 0.00002 -num 1000 -S 4 4"<<endl;
	//cout<<"hybrid-Lambda -spng 2_tax_sp1_gener_num -pop_size 10000 -mu 0.0002 -mm 1.4 -seg -num 100 -S 4 4"<<endl;
	//cout<<endl;
//}
