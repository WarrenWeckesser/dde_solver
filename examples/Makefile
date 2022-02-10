# This is a GNU Makefile.

FC = gfortran
DDE_SOLVER_DIR = ..
PROGS = b2 b2g c2g ex4p4p1 ex4p4p2 ex4p4p3 ex4p4p4 ex4p4p5 ex4p7 ex4p8 exsd1 \
        infection laserex neutvan secdelay volinteq

all: dde_solver_m.o $(PROGS)

dde_solver_m.o: $(DDE_SOLVER_DIR)/dde_solver_m.f90
	$(FC) -c $(FFLAGS) $<

b2: b2.f90 dde_solver_m.o
	$(FC) $^ -o $@

b2g: b2g.f90 dde_solver_m.o
	$(FC) $^ -o $@

c2g: c2g.f90 dde_solver_m.o
	$(FC) $^ -o $@

ex4p4p1: ex4p4p1.f90 dde_solver_m.o
	$(FC) $^ -o $@

ex4p4p2: ex4p4p2.f90 dde_solver_m.o
	$(FC) $^ -o $@

ex4p4p3: ex4p4p3.f90 dde_solver_m.o
	$(FC) $^ -o $@

ex4p4p4: ex4p4p4.f90 dde_solver_m.o
	$(FC) $^ -o $@

ex4p4p5: ex4p4p5.f90 dde_solver_m.o
	$(FC) $^ -o $@

ex4p7: ex4p7.f90 dde_solver_m.o
	$(FC) $^ -o $@

ex4p8: ex4p8.f90 dde_solver_m.o
	$(FC) $^ -o $@

exsd1: exsd1.f90 dde_solver_m.o
	$(FC) $^ -o $@

infection: infection.f90 dde_solver_m.o
	$(FC) $^ -o $@

laserex: laserex.f90 dde_solver_m.o
	$(FC) $^ -o $@

neutvan: neutvan.f90 dde_solver_m.o
	$(FC) $^ -o $@

secdelay: secdelay.f90 dde_solver_m.o
	$(FC) $^ -o $@

q1damodule.o: q1damodule.f90
	$(FC) -c $(FFLAGS) $< -o $@
volinteq: volinteq.f90 q1damodule.o dde_solver_m.o
	$(FC) $^ -o $@

clean:
	rm -f $(PROGS) *.o *.mod
	rm -f *.dat debug.ans
