gfortran -c ../../dde_solver_m.f90
gfortran -c ../testing.f90
gfortran -c feldstein_neves1984.f90
gfortran feldstein_neves1984_test.f90 feldstein_neves1984.o testing.o dde_solver_m.o -o feldstein_neves1984_test
./feldstein_neves1984_test
