SRCPATH   = ../src/
LIBPATH   = ../lib/
LIB       = optim
STATUT=debug
OPTS      =  -Wall -fPIC 
ifeq ($(STATUS), debug)
OPTS     += -g -DDEBUG
else
OPTS     += -O2
endif
IPLP = $(SRCPATH)affinescaling.o $(SRCPATH)linprosimpx.o $(SRCPATH)initlp.o \
       $(SRCPATH)linprosimp.o $(SRCPATH)solvetriglp.o \
       $(SRCPATH)modchol.o  $(SRCPATH)modchol2.o  
$(IPLP):$(FRC)

      

FRC:
	@FRC=$(FRC)	

all: EXE

$(LIB): $(IPLP)
	ar cr $(LIBPATH)$@.a $(IPLP)
	ranlib $(LIBPATH)$@.a
	
	
EXE: $(LIB)
	gfortran $(OPTS) -o ../bin/linprog  $(SRCPATH)main.F90 ../lib/optim.a  -llapack
	gfortran $(OPTS) -o ../bin/modchol  $(SRCPATH)main_modchol.F90 ../lib/optim.a  -llapack

clean:
	rm -rf  ../src/*.o
	
cleanall: clean 
	rm -rf  $(LIB).a $(LIB).so ../bin/linprog
	
%.o: %.F90
	gfortran $(OPTS) -c $< -o $@		
