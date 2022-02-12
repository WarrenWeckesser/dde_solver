dde_solver
==========

This repository contains an updated version of the DDE_SOLVER Fortran
code written by S. Thompson and L. Shampine [1].

Major changes to dde_solver_m.f90
---------------------------------
* *Use assumed-shape arrays in all interfaces.*
  This fixed the segmentation faults that were occurring in the old version.
  (Users' old code will have to be updated to match the new interfaces.)
* *Removed ODEAVG,* and removed all the code in dde_solver_m.f90 related to
  averaging.  (This was done on the recommendation of S. Thompson.)
* *Added the function DDE_SOLVER_VERSION* that returns the software version.

The `examples/` directory has S. Thompson's original examples (with some
minor updates).

New unit tests are maintained in the `tests/` directory.

-----

1.  S. Thompson, L.F. Shampine, A friendly Fortran DDE solver,
    *Applied Numerical Mathematics*, Volume 56, Issues 3–4, 2006,
    pages 503-516, ISSN 0168-9274,
    https://doi.org/10.1016/j.apnum.2005.04.027.
