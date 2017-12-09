dde_solver
==========

This repository contains my work on updating the DDE_SOLVER Fortran
code originally written by L. Shampine and S. Thompson:

    http://www.radford.edu/~thompson/ffddes/

When the DDE_SOLVER code is compiled with recent versions of gfortran,
running the solver may result in segmentation faults or other crashes.
The problem appears to be that not enough interface declarations have
been made in the code, so gfortran does not know that some array arguments
are assumed-shape.  It passes invalid data in some subroutine calls,
resulting in a segmentation fault.

The initial commit adds the file "dde_solver_m.f90".  This file was originally
called "dde_solver_m_unix.f90" in the ZIP file available on S. Thompson's web
page.  The only difference between the original "dde_solver_m.f90" and
"dde_solver_m_unix.f90" was the line endings.  For this repository, I am
maintaining just one file, with "Unix" (i.e. practically every system
except Windows) line endings.

Major changes
-------------
* Use assumed-shape arrays in all interfaces.  (Old code will have to be updated to
  match the new interfaces.)
* Removed ODEAVG, and removed all the code in dde_solver_m.f90 related to averaging.
* Added the function DDE_SOLVER_VERSION that returns the software version.
