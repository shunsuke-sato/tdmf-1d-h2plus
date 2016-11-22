#FC = gfortran -O2 ## gfotran
FC = mpiifort -O3 -ipo -ip -xHOST ##cal*@nucl_thr w/o openMP

#LN = -llapack -lblas #other
LN = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm



VPATH = src:object
SRC = $(shell cd src ;ls *.f90 ;cd ..)
OBJ = $(SRC:.f90=.o)
OBJ_dir = $(addprefix object/,$(OBJ))

PROG = tdmf

$(PROG):global_variables.o $(OBJ)
	$(FC) -o $(PROG) $(OBJ_dir) $(LN)

main.o:main.f90
	$(FC) -c $< $(LN);mv $@  object 
%.o:%.f90
	$(FC) -c $< $(LN);mv $@  object 


clean:
	rm  -f  object/*.o  *.mod tdmf
clean_complete:
	rm  -f *~  */*~ */*/*~ object/*.o  */*.mod *.mod tdmf */#* *.out log.log
