#
# This file is part of ssparse.
#
# Copyright (C) 2012, Technical University of Denmark
# https://github.com/nefan/ssparse
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

# C flags
CF += -std=c99 

# avx? sse4?
CF += -mavx
# CF += -msse4

# using CUDA?
# USE_CUDA = true

ifdef USE_CUDA

CF := $(CF) -DSSP_CUDA

# CUDA compiler
NVCC = nvcc

NVCF  = $(CF)
NVCF := $(shell echo $(NVCF) | sed "s/-std=c99 //g")
NVCF := $(shell echo $(NVCF) | sed "s/-mavx //g")
NVCF := $(shell echo $(NVCF) | sed "s/-fexceptions //g")
NVCF := $(shell echo $(NVCF) | sed "s/-fPIC //g")
NVCF := $(NVCF) -Xcompiler -fpic

endif
