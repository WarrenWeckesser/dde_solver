gfortran -c ../../dde_solver_m.f90
gfortran -c ../testing.f90
gfortran -c paul1992.f90
gfortran paul1992_test.f90 paul1992.o testing.o dde_solver_m.o -o paul1992_test
./paul1992_test
