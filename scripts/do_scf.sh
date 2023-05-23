for i in `seq 1 3000` ; do
echo "----------" $i "-----------"
  if [ -f "negf.converged" ] || [ -f negf.stop ] ; then
    exit
  else
    ./scf.sh > scf.out
    err=`echo $?`
    if [ $err -ne 0 ] ; then
      echo "SCF failed " $err
      exit $err
    fi

  fi
done
