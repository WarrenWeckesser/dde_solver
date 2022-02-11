gfortran -g -c ../../dde_solver_m.f90
gfortran -g -c ../testing.f90
gfortran -g -c delay_oscillator.f90
gfortran -g delay_oscillator_test.f90 delay_oscillator.o testing.o dde_solver_m.o -o delay_oscillator_test
./delay_oscillator_test
