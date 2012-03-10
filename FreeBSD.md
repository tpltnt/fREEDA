Notes for FreeBSD
=================

Ports used
----------
- [g95](http://www.freshports.org/lang/g95/) - Fortran compiler: cd /usr/ports/lang/g95/ && make install clean

  - does not build with FreeBSD clang version 3.0 (branches/release_30 142614) 20111021
  
- [BLAS](http://www.freshports.org/math/blas/) - Basic Linear Algebra Subroutines: cd /usr/ports/math/blas/ && make install clean
- [LAPACK](http://www.freshports.org/math/lapack/) - A library of Fortran 77 subroutines for linear algebra: cd /usr/ports/math/lapack/ && make install clean

useful links
------------

- http://www.netlib.org/lapack95/
- http://gcc.gnu.org/fortran/
- http://gcc.gnu.org/wiki/GFortran
- http://gfortran.com/

status
------

- not in ports collection
- building gfortran: [wiki instructions](http://gcc.gnu.org/wiki/GFortranSource) and [Installing GCC](http://gcc.gnu.org/install/)

  1. download source
  ```wget http://www.artfiles.org/gnu.org/gcc/gcc-4.2.1/gcc-fortran-4.2.1.tar.bz2```

  2. unpack it
  ```tar -xf gcc-fortran-4.2.1.tar.bz2```

  3. change into source directory
  ```cd gcc-4.2.1```

- building trilinos fails

  - /usr/local/bin/g95 CMakeFiles/cmTryCompileExec.dir/testFortranCompiler.f.o -o cmTryCompileExec
  - ld: cannot find -lf95
  - gmake[1]: *** [cmTryCompileExec] Error 1
  - gmake: *** [cmTryCompileExec/fast] Error 2
