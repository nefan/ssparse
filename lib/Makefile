#
# This file is part of ssparse.
#
# Copyright (c) 2006, Timothy A. Davis.
# File derived from:
# CXSparse: a Concise Sparse matrix package - Extended.
# http://www.suitesparse.com
# 
# Modified by Stefan Sommer (shso@elektro.dtu.dk), November 2012
#
# ssparse is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# ssparse is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this Module; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# 

default: all

include ../../SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk
include ../ssparse_config.mk

I = -I../../SuiteSparse/SuiteSparse_config -I ../../SuiteSparse/CXSparse/Include -I../include -I$(INSTALL_INCLUDE)

all: libssparse.a

SSP_SOURCE = ssp_splsolve.c ssp_sputsolve.c ssp_spdot.c

SSP_DI_OBJ =  ssp_splsolve_di.o ssp_sputsolve_di.o ssp_spdot_di.o

SSP_DL_OBJ =  ssp_splsolve_dl.o ssp_sputsolve_dl.o ssp_spdot_dl.o

SSP_CI_OBJ =  ssp_splsolve_ci.o ssp_sputsolve_ci.o ssp_spdot_ci.o

SSP_CL_OBJ =  ssp_splsolve_cl.o ssp_sputsolve_cl.o ssp_spdot_cl.o

SSP_CUDA = ssp_cuda.o ssp_spdot_cuda.o

SSP = $(SSP_DI_OBJ) $(SSP_DL_OBJ) $(SSP_CI_OBJ) $(SSP_CL_OBJ) 
ifdef USE_CUDA
	SSP+= $(SSP_CUDA)
endif

$(SSP): ../include/ssparse.h Makefile

%_di.o : ../src/%.c
	$(CC) $(CF) $(I) -c $< -o $@

%_dl.o : ../src/%.c
	$(CC) $(CF) $(I) -DCS_LONG -c $< -o $@

%_ci.o : ../src/%.c
	$(CC) $(CF) $(I) -DCS_COMPLEX -c $< -o $@

%_cl.o : ../src/%.c
	$(CC) $(CF) $(I) -DCS_LONG -DCS_COMPLEX -c $< -o $@

%_cl.s : ../src/%.c
	$(CC) $(CF) $(I) -DCS_LONG -DCS_COMPLEX -S -fverbose-asm $< -o $@

%_cuda.o : ../src/%_cuda.cpp
	# $(NVCC) $(NVCF) $(I) -DCS_COMPLEX -c $< -o $@
	$(CC) $(CF) $(I) -DCS_LONG -DCS_COMPLEX -c $< -o $@

libssparse.a: $(SSP)
	$(ARCHIVE)  libssparse.a $(SSP)
	- $(RANLIB) libssparse.a

clean:
	- $(RM) *.o
	- $(RM) *.s

purge: distclean

distclean: clean
	- $(RM) *.a
