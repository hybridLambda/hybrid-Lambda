
Hybrid-Lambda is a software package that can simulate gene trees within a rooted
species network or a rooted species tree under the coalescent process. The
main feature of this program is that users can choose to use the standard
Kingman coalescent process, which produces bifurcating genealogies, or two
other Lambda coalescent processes, which produce multifurcating genealo-
gies. The other feature is that hybrid sim uses extended Newick formatted
strings to make it easier to represent hybridization events between species.


##INSTALLATION
To install hybrid-Lambda, first install the following packages and libraries

on Debian/Ubuntu based systems:
```bash
apt-get install build-essential autoconf autoconf-archive doxygen graphviz
```
on Mac OS:
```bash
port install automake autoconf autoconf-archive doxygen graphviz
```

then type the following commands:
```bash
./bootstrap
make
```


##ASSUMPTION
Input network files are written in extended newick format.
 

##LICENCE
You can freely use all code in this project under the conditions of the GNU
GPL Version 3 or later.

##HOW IT WORKS

Program parameters and options:

Options              | Useage
-------------------- | -----------------------
-h or -help          | Help. List the following content.
-spcu STR          | Input the species network/tree string through command line or a file. Branch lengths of the INPUT are in coalescent unit.
-spng STR          | Input the species network/tree string through command line or a file. Branch lengths of the INPUT are in number of generation.
-pop STR/FLT           | Population sizes are defined by a single numerical constant, or a string which specifies the population size on each branch. The string can be input through command line or a file. By **default**, population size 10,000 is used.
-mm STR/FLT            | Multiple merger parameters are defined by a single numerical constant, or a string which speifies the parameter on each branch. The string can be input through command line or a file. By **default**, Kingman coalescent is used.
-S INT INT ...         | Specify the number of samples for each taxon.
-num INT               | The number of gene trees will be simulated.
-seed INT           | User define random SEED
-mu FLT               | User defined constant mutation rate per locus. By **default** mutation rate 0.00005 is used.
-o STR [option]   | Specify the file name prefix for simulated gene trees. Prefix is set as "OUT" by **default**, When options are not specified, only output trees with branch lengths are in coalescent unit.
-sim_mut_unit        | Convert the simulated gene tree branch lengths to mutation unit.
-sim_num_gener       | Convert the simulated gene tree branch lengths to number of generations.
-sim_num_mut         | Simulate numbers of mutations on each branch of simulated gene trees.
-sim_Si_num          | Generate a table, which includes the number of segregating sites and the total branch length of the gene tree, as well as the TMRCA.
-f                   | Generate a topology frequency table of a set of input trees or simulated gene trees. 
-gt STR             | Specify the FILE NAME of trees to analyse tree topology frequencies.
-seg  |  Generate segregating site data.
-mt STR  |  Specify the FILE NAME of trees to generate segregating site data. Tree branch lengths indicate number of mutations on the branch.
-mono                | Generate a frequency table of monophyletic, paraphyletic and polyphyletic trees. 
-plot/-dot [option]  | Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree).
-branch              | Branch lengths will be labelled in the figure.
-plotF/-dotF FILE    | Generated figure will be saved in FILE.

##Examples:
```bash
hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 3 -seed 2 -o example1
hybrid-Lambda -spcu trees/4_tax_sp_nt1_para -o example2 -num 2 -mu 0.00003 -sim mut unit -sim num mut
hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num N -pop 25000 -sim num gener
hybrid-Lambda -spng '(A:50000,B:50000)r;' -pop '(A:50000,B:50000)r:40000;'
hybrid-Lambda -spcu '((((A:1.1,B:1.1):2.1,a:2.2):1.1,13D:.2):.3,4:.3);' -S 2 4 3 6 5
hybrid-Lambda -spcu '(A:1,B:1)r;' -mm '(A:1.9,B:.2)r:2;' -S 3 4
hybrid-Lambda -spcu trees/7_tax_sp_nt1_para -dot -branch
hybrid-Lambda -spcu trees/4_tax_sp1 -num 1000 -o GENE_TREE_FILE -f
hybrid-Lambda -spcu trees/4_tax_sp1 -num 1000 -o GENE_TREE_FILE -fF FRENQUENCY_FILE
hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 1000 -o GENE -fF OUTPUT
hybrid-Lambda -gt GENE_coal_unit -f 
hybrid-Lambda -spcu '(A:5,B:5)r;'-mono -num 100 -mm .1 -S 4 4
```
