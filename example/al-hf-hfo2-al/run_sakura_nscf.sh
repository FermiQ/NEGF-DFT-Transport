#$ -pe q*  32
#$ -cwd

module load gcc/gcc-9.2.1 intel/intel-latest intelmpi/intelmpi-latest


export PETSC_DIR=$HOME/prog/petsc/petsc3.14_intelmpi_opt
export LD_LIBRARY_PATH=$PETSC_DIR/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

export ppn=32

export I_MPI_EXTRA_FILESYSTEM=1


export solver="superlu_dist"
export transomat_bin=$HOME/prog/transomat-release/bin/transomat
export conquest_bin=$HOME/prog/conquest/CONQUEST-release/bin/Conquest

cd ecc

mv Conquest_input Conquest_input.dft
cat Conquest_input.dft  | sed 's/negf.LoadKnegf f/negf.LoadKnegf t/g' > Conquest_input.tmp

../scripts/do_scf.sh > do_scf.out 2>do_scf.err

cd ..
