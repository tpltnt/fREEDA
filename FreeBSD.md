Notes for FreeBSD
=================

Ports used
----------
- [g95](http://www.freshports.org/lang/g95/) - Fortran compiler: cd /usr/ports/lang/g95/ && make install clean

  - does not build with FreeBSD clang version 3.0 (branches/release_30 142614) 20111021
  
- [BLAS](http://www.freshports.org/math/blas/) - Basic Linear Algebra Subroutines: cd /usr/ports/math/blas/ && make install clean
- [LAPACK](http://www.freshports.org/math/lapack/) - A library of Fortran 77 subroutines for linear algebra: cd /usr/ports/math/lapack/ && make install clean
