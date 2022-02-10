gfortran -c ../../dde_solver_m.f90
gfortran -c ../testing.f90
gfortran -c delay_oscillator.f90
gfortran delay_oscillator_test.f90 delay_oscillator.o testing.o dde_solver_m.o -o delay_oscillator_test
./delay_oscillator_test
