# This is a GNU Makefile.

FC = gfortran
DDE_SOLVER_DIR = ..
SINGLE_FILE_PROGS = b2 b2g c2g ex4p4p1 ex4p4p2 ex4p4p3 ex4p4p4 ex4p4p5 ex4p7 \
                    ex4p8 exsd1 infection laserex neutvan secdelay
ALL_PROGS = $(SINGLE_FILE_PROGS) volinteq

all: dde_solver_m.o $(ALL_PROGS)

dde_solver_m.o: $(DDE_SOLVER_DIR)/dde_solver_m.f90
	$(FC) -c $(FFLAGS) $<

$(SINGLE_FILE_PROGS): %: %.f90 dde_solver_m.o
	$(FC) $(FFLAGS) $^ -o $@

q1damodule.o: q1damodule.f90
	$(FC) -c $(FFLAGS) $< -o $@
volinteq: volinteq.f90 q1damodule.o dde_solver_m.o
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm -f $(ALL_PROGS) *.o *.mod
	rm -f *.dat debug.ans
