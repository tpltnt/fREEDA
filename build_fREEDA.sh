#!/bin/bash
# Run this script from the top level freeda-<version> directory

echo "****************************************************"
echo "*"							 
echo "*" Welcome to the fREEDA installation script        
echo "*"							 
echo "****************************************************"

# First make the libraries
cd libs

echo "****************************************************"
echo "*"						
echo "*" The FFTW library   
echo "*"					
echo "****************************************************"

cd fftw-2.1.5
./configure
make 
if [ "$?" = 2 ]; then
  echo "***FFTW compilation failed"
  exit 1
else  
  echo "FFTW library finished"
  cp fftw/.libs/libfftw.a ../../simulator/lib
  cp rfftw/.libs/librfftw.a ../../simulator/lib
  cd ../../
fi  


if [ "$?" = 0 ]; then
  echo "****************************************************"
  echo "*"						
  echo "*" "fREEDA support libraries compilation COMPLETE!"
  echo "*" "Now beginning compilation of fREEDA routines"
  echo "*"					
  echo "****************************************************"
fi

#export my_path=`pwd`

# Insert correct path in simulator/makefile.defs
#sed s+freeda_dir_path+$my_path+ simulator/makefile.defs > simulator/makefile.defs.tmp
#if [ -f simulator/makefile.defs.tmp ]; then
#  mv simulator/makefile.defs.tmp simulator/makefile.defs
#  rm simulator/makefile.defs.tmp
#fi


cd simulator

make dep

if [ "$?" = 0 ]; then
  make
  if [ "$?" = 0 ]; then
    echo "fREEDA COMPILATION COMPLETE"
  else
    echo "***Make  error."
    exit 1
  fi
else 
  echo "Dependency error!"
  echo "fREEDA compilation failed!"
  exit 1
fi


