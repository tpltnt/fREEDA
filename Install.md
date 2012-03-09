Last updated Nov 2010.

This file contains information related to installation of fREEDA and
ifreeda which is a GUI front-end to fREEDA. The README and
README _ GUI files contain useful usage information.

The use of the '$' symbol denotes the shell prompt.

System Requirements
-------------------
- GNU Make
- GNU C compiler gcc/g++ version 4.4.x 
- GNU fortran compiler gfortran 4.4.x
- GNU flex-2.5.xx and bison-2.4.x
- The bash shell version 4.1.5
- CMake 2.8.x for the Trilinos installation

For ifreeda

- Qt 4.5.0. This includes the qmake tool, and the associated toolkit library.

It is encouraged to use Qt 4.5.0 only. Newer versions of Qt have not been
tested.

Directory Structure
-------------------
freeda-<version-number> (main directory)

	- libs (containes the libraries that are linked with freeda)
	- doc (freeda PDF documentation)
	- ifreeda (contains GUI code)
	- simulator (the main freeda code)
	
		- analysis 
		- compat
		- elements
		- inout
		- lib
		- network
		- projects (contains examples to test with the GUI)

	- test (contains some sample netlists)


Compiling and Installing fREEDA
-------------------------------
The first step is installing the Trilinos libraries and associated
dependencies. Since we do not maintain Trilinos, we cannot provide
support for installing it. However, we will provide the steps we took to
installing Trilinos. Look at the official Trilinos instructions if you
want to use options in addition to what we used. In order to get fREEDA
working, you have to at least compile the portions of Trilinos that we use.

TRILINOS INSTALLATION

First ensure that BLAS (http://www.netlib.org/blas/) and Lapack
(http://www.netlib.org/lapack/) are installed, make sure to compile them
with the gfortran compiler. For AMD64 platforms, use the -fPIC option when
compiling BLAS and Lapack.

For BLAS, unpack the blas.tgz archive and edit the file make.inc.
Change all the g77 instances to gfortran. The run make, and cp the generated
library file to /usr/lib.

```$ make```

Copy the library file to the system-level folder. This makes it easy for
Trilinos to find it.
```$ sudo cp blas_LINUX.a /usr/lib/libblas.a```

For LAPACK, unpack the lapack.tgz archive and use the commands
```$ cp make.inc.example make.inc
$ make lapacklib```

Copy the library file to the system-level folder. This makes it easy for
Trilinos to find it.
```$ sudo cp lapack_LINUX.a /usr/lib/liblapack.a```

Ensure that g77 is not installed on your system. This will cause problem
with the Trilinos installation.

Download Trilinos version 10.0.2 from http://trilinos.sandia.gov/download/trilinos-10.0.html.
There are newer versions of Trilinos available but we have not tested fREEDA
with them yet.

Unpack the Trilinos archive into the /opt directory. This will create a
directory called /opt/trilinos-10.0.2-Source. Note that for this and subsequent
operations, you may need admin rights on your computer.

Trilinos is designed to do both serial and parallel computations. To install
the serial capabilities of Trilinos, create a directory under trilinos-10.0.2-Source
called LINUX_SERIAL. The directory structure now looks like this:
/opt/trilinos-10.0.2-Source/LINUX_SERIAL.
Inside the LINUX_SERIAL directory create a file named do-configure and it
should contain the following lines:

```
#!/bin/sh
# cmake script for a serial build of trilinos on 64-bit Linux

export TRILINOS_HOME=/opt/trilinos-10.0.2-Source

EXTRA_ARGS=$@

cmake \
  -D CMAKE_BUILD_TYPE:STRING=DEBUG \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
  -D Trilinos_ENABLE_Komplex:BOOL=ON \
  -D Trilinos_ENABLE_Sacado:BOOL=ON \
  -D Trilinos_ENABLE_Triutils:BOOL=ON \
  -D CMAKE_INSTALL_PREFIX:PATH=$TRILINOS_HOME/LINUX_SERIAL \
  $EXTRA_ARGS \
  ${TRILINOS_HOME}
```


This file contains all the Trilinos packages required by fREEDA. Save this
file and make it executable by issuing the command
```$ chmod +x do-configure```

Run this executable file.
```$ ./do-configure```

Then compile Trilinos.
```$ make```

When the Trilinos compilation  is finished, go to your freeda-<version-number>/simulator
directory.

Open the file makefile.defs in a text editor and update the variables FRHOME,
TSRC and TLIB to indicate the path of your fREEDA and Trilinos installations
respectively.
FRHOME is the path of your fREEDA home directory. So if you unpacked the fREEDA
tarball into your home directory, this variable would be
FRHOME = /home/myusername/freeda-2.0 (for version 2.0 of fREEDA)
For example, if you have followed the instructions above, these variables would look like
TSRC = /opt/trilinos-10.0.2-Source/packages
TLIB = /opt/trilinos-10.0.2-Source/LINUX_SERIAL/packages

Save and exit this file, and go one level up to the freeda-<version-number> directory.
Run the script build_fREEDA.
``$ cd ..
$ ./build_fREEDA```

This will compile all of the fREEDA code. Once this completes succesfully,
the freeda binary will be created in freeda-<version-number>/simulator. From
there, you can run the simulator on any of the sample netlists in the test directory.
For example, from freeda-<version-number>/simulator:
```./freeda <netlist_name>```


Compiling and Installing the GUI (ifreeda)
------------------------------------------

Acknowledgements: 
The GUI packaged with fREEDA is based on a gui used in the software simulator Qucs. The Qucs project can be found at http://qucs.sourceforge.net/
The GUI is based on the Qt toolkit (http://qt.nokia.com)

Before you can compile the GUI packaged with fREEDA you must do a few things first.

1. Get Qt 4.5.0
Download and install Qt version 4.5.0. You can find a package suitable for your platform at:
ftp://ftp.trolltech.com/qt/source/
For example, on our Linux installation, we used the file qt-x11-opensource-src-4.5.0.tar.gz. 
This is the only version that has been tested with the GUI.

2. Install Qt 4.5.0
Follow the INSTALL instructions for Qt. Specifically setting up the environment variables for Qt. During the configure phase, we used the following command:
./configure -static -fast
This ensures that a static Qt library file is built.

The fREEDA GUI  will need to know the location of the Qt includes, libraries and qmake. These paths need to be in your PATH. The instructions for doing this properly are in the Qt INSTALL documentation.

3. Compiling ifREEDA
In the freeda-<version-number>/ifreeda directory, you must edit a file called ifreeda.pro. In that file are two variables called TSRC and TLIB. These point to your Trilinos library installation. 
In our case, we have Trilinos 10.0.2 installed in the directory /opt. So these variables are named as
TSRC = /opt/trilinos-10.0.2-Source/packages
TLIB = /opt/trilinos-10.0.2-Source/LINUX_SERIAL/packages
You might choose a different path, in which case you will update these variables accordingly.

Run the ./build_GUI script, which you will find in the freeda-<version-number> directory.

This will create the ifreeda binary in the freeda-<version-number>/ifreeda directory.
To run this, from the ifreeda directory, type
$ ./ifreeda


Contact information
-------------------
This INSTALL file and the compilation process is maintained by Nikhil Kriplani.
Please direct any issues and queries to nkriplani-HERE-COMES-THE-AT-SIGN-ncsu.edu. Getting the system to compile correctly is always a work in progress, so any feedback to help improve 
the process is welcome. If the compilation is succesful, send a friendly note and let us know.
