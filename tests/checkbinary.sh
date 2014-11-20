#!/bin/bash
function test_hybrid-Lambda {
  echo -n " hybrid-Lambda $@ "
  for i in `seq 1 1`; do
    echo -n "."

    # Test using hybrid-Lambda self-checks
    ./hybrid-Lambda_dbg $@ -seed $i > /dev/null 
    if [ $? -ne 0 ]; then
      echo ""
      echo "Executing \"./hybrid-Lambda_dbg $@ -seed $i\" failed."
      echo "Debug Call: make -mj2 hybrid-Lambda_dbg && ./hybrid-Lambda_dbg $@ -seed $i 2>&1 | less"
      exit 1
    fi

    # Test for memory leaks
    valgrind --error-exitcode=1 --leak-check=full -q ./hybrid-Lambda $@ -seed $i > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Valgrind check of \"./hybrid-Lambda $@ -seed $i\" failed."
      exit 1
    fi

  done
  echo " done."
}

echo "Testing Examples"
	test_hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -print
	test_hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 3 -seed 2 -o example1 || exit 1
	test_hybrid-Lambda -spcu trees/4_tax_sp_nt1_para -o example2 -num 2 -mu 0.00003 -sim_mut_unit -sim_num_mut || exit 1
	test_hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 100 -pop 25000 -sim_num_gener || exit 1
	test_hybrid-Lambda -spng '(A:50000,B:50000)r;' -pop '(A:50000,B:50000)r:40000;' || exit 1
	test_hybrid-Lambda -spcu '((((A:1.1,B:1.1):2.1,a:2.2):1.1,13D:.2):.3,4:.3);' -S 2 4 3 6 5 || exit 1
	test_hybrid-Lambda -spcu '(A:1,B:1)r;' -mm '(A:1.9,B:.2)r:2;' -S 3 4 || exit 1
	test_hybrid-Lambda -spcu trees/7_tax_sp_nt1_para -dot -branch || exit 1
	test_hybrid-Lambda -spcu trees/4_tax_sp1.tre -num 1000 -o GENE -f || exit 1
	test_hybrid-Lambda -gt GENE_coal_unit -f  || exit 1
    test_hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -o GENE -num 100 -mu 0.00003 -sim_num_mut || exit 1
	test_hybrid-Lambda -mt GENE_num_mut -seg  || exit 1
	test_hybrid-Lambda -spcu '(A:5,B:5)r;' -mono -num 100 -mm .1 -S 4 4 || exit 1
echo ""
