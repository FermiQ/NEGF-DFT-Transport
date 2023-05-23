.SUFFIXES: .f90

FC = ifort

FFLAGS = -g -traceback -c

LINK = -L.




OBJS = make_conquest_mod.o transomat_make_conquest_input.o
%.o %.mod : %.f90
	$(FC) $(FFLAGS)  $*.f90

transomat_make_conquest_input: $(OBJS)
	$(FC)  -o $@ $(OBJS) $(LINK)


clean:
	rm -f *.o *.mod transomat_make_conquest_input

