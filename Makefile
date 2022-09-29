LIB_LIST = -L${BLASDIR} ${BLASLIB} -lpthread

F95 = gfortran
OPT = -O3 -Wall

all: grayscott

create_matrices.o: create_matrices.f90
	$(F95) -c $(OPT) $<

create_banded.o: create_banded.f90
	$(F95) -c $(OPT) $<

banded_crank.o: banded_crank.f90
	$(F95) -c $(OPT) $<

test_banded.o: test_banded.f90
	$(F95) -c $(OPT) $<

heun.o: heun.f90
	$(F95) -c $(OPT) $<

save_fields.o: save_fields.f90
	$(F95) -c $(OPT) $<

initial.o: initial.f90
	$(F95) -c $(OPT) $<

timestepping_improved.o: timestepping_improved.f90
	$(F95) -c $(OPT) $<

grayscott.o: grayscott.f90
	$(F95) -c $(OPT) $<

grayscott: grayscott.o timestepping_improved.o initial.o save_fields.o heun.o test_banded.o banded_crank.o create_matrices.o create_banded.o
	$(F95) -o grayscott $^ $(LIB_LIST)

clean:
	rm -f *.o core.* grayscott
