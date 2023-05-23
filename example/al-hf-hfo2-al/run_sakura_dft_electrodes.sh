#$ -pe q*  32
#$ -cwd

module load gcc/gcc-9.2.1 intel/intel-latest intelmpi/intelmpi-latest

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

export ppn=12

export I_MPI_EXTRA_FILESYSTEM=1

export conquest_bin=$HOME/prog/conquest/CONQUEST-release/bin/Conquest

cd left_electrode

mpiexec.hydra -ppn $ppn -genvall -envall -launcher ssh \
$conquest_bin > conquest.out 2>conquest.err

cd ..
cd right_electrode

mpiexec.hydra -ppn $ppn -genvall -envall -launcher ssh \
$conquest_bin > conquest.out 2>conquest.err

cd ..
