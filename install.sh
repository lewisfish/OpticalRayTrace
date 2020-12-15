#!/bin/bash

function showhelp
{

  echo 'Usage: ./install.sh [option] [option]...'
  echo 'Compile and run mcgrid.'
  echo 
  echo '   -h, --help            Shows this dialog.'
  echo '   -n, --threads         Compiles and run code on n threads. Default is 32 threads.'
  echo '   -d, --debug           Compiles and runs code with debug flags on n threads.'
  echo '   -m, --make            Compiles the code with warning enabled.'

}

function makebuild
{
  if [ "$comp" = 'gnu' ];then
      string="FCOMP=gfortran"
  fi

  if [ "$debug" = 1 ];then
    make clean && make debug $string
  elif [ "$make" = 1 ];then
    make clean && make build $string
  else
    if [ "$NUM_THREADS" = 1 ];then
      make clean && make $string
    else
      make clean && make mp $string
    fi
  fi
}

function createdirs
{

  if [ ! -d "build" ]; then
     mkdir "build"
  fi
  cd build
  ndirec="$(pwd)"
  cd ..
  if [ ! -d "bin" ]; then
     mkdir "bin"
  fi
  cd bin
  bdirc="$(pwd)"
  cd ..
  cd src
}

function run
{
  for i in *; do
     if [ "${i}" != "${i%.mod}" ];then
        cp "${i}" "$ndirec"
     fi
     if [ "${i}" != "${i%.o}" ];then
        mv "${i}" "$ndirec"
     fi
  done


  if [ "$make" = "1" ]; then #just make code
      exit 0
  fi

  mv raytrace "$bdirc" && echo " "&& echo "*****Install complete*****" && echo " "

  clear
  cd ../bin

  ./raytrace
}

#defaults
NUM_THREADS=32
debug=0
help=0
comp="gnu"
make=0

set -e

createdirs

while [ "$1" != "" ]; do
    case $1 in
        -n | --threads )        NUM_THREADS=$2
                                ;;
        -h | --help )           showhelp
                                exit
                                ;;
        -m | --make )           make=1
                                makebuild
                                exit
                                ;;
        -d | --debug )          debug=1
                                ;;
    esac
    shift
done

makebuild
if [ "$NUM_THREADS" > 1 ]; then
  export OMP_NUM_THREADS=$NUM_THREADS
fi
run