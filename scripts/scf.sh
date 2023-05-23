
rm  subdiv* -f

mpiexec.hydra -ppn $ppn -genvall -envall -launcher ssh $transomat_bin \
-log_view :log \
-pc_factor_mat_solver_type $solver
-mat_superlu_dist_fact SamePattern_SameRowPerm -mat_superlu_dist_rowperm NOROWPERM \
-mat_superlu_dist_colperm METIS_AT_PLUS_A > trans.out 2> trans.err

err=`echo $?`
if [ $err -ne 0 ] ; then
  echo "transomat failed " $err
  exit 1
fi

if [ ! -f negf.converged ] ; then

  mpiexec.hydra -ppn $ppn -genvall -envall -launcher ssh \
  $conquest_bin > conquest.out 2>conquest.err

  err=`echo $?`
  if [ $err -ne 0 ] ; then
    echo "Conquest failed " $err
    exit 2
  fi


else
 echo "NEGF converged"
fi
