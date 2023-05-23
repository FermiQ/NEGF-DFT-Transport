transomat_root=/home/prog/ownstuff/fortran/programs/transomat_dev/release_v0.1/transomat-release

$transomat_root/tools/transomat_make_conquest_input/transomat_make_conquest_input coord.in 14.0 110.0

scp pseudo/Al.ion left_electrode
scp pseudo/Al.ion right_electrode
scp pseudo/*.ion ecc

scp transp.ini ecc
