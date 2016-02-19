FC = gfortran
FCFLAGS = -g -fbounds-check -ffixed-line-length-132
FCFLAGS = -O2
# FFLAGS= -Mextend -O3 -Munroll -Mnoframe # used this for normal stuff
 FFLAGS = -g -fbounds-check -ffixed-line-length-132
# FFLAGS = -g  -132 -kind=byte
# flags forall (e.g. look for system .mod files, required in gfortran)
# FCFLAGS += -I/usr/include
# LIBS = libmyblas.a libmylapack.a libmyblas.a 
all: ujorxn

# links
#logSumDiff.o: mrgrnk.mod refsor.mod
main.o: rxn.o 
multLogMat.o: logSumDiff.o paramsrxn.o l2d.o
rxn.o: multLogMat.o logSumDiff.o lpopul.o paramsrxn.o l2d.o
logSumDiff.o: mrgrnk.o refsor.o 
l2d.o: logSumDiff.o

# general rules
ujorxn: main.o rxn.o logSumDiff.o lpopul.o multLogMat.o paramsrxn.o mrgrnk.o refsor.o l2d.o
	$(FC) $(FCFLAGS) -o  $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

%.o: %.F90
	$(FC) -c $<

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)

