#!/bin/bash

export FC=mpiifort
export CC=mpiicc
export F77=mpiifort

export FCFLAGS="`nc-config --flibs` -I`nc-config --includedir` -I/home/Mikyung.Lee/am5phys_cleanup/FMS/fms_intel/include -O0 -g -trackback -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -i4 -nowarn -L/home/Mikyung.Lee/am5phys_cleanup/FMS/fms_intel/lib"
export LDFLAGS="-L/home/Mikyung.Lee/am5phys_cleanup/FMS/fms_intel/lib"
export CFLAGS="-O0 -sox -qopenmp -debug minimal `nc-config --cflags`"
