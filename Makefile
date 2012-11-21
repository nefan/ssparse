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

#------------------------------------------------------------------------------
# ssparse Makefile
#------------------------------------------------------------------------------

VERSION = 0.0.1

include ssparse_config.mk

default: C

include ../SuiteSparse_config/SuiteSparse_config.mk

C:
	( cd lib ; $(MAKE) )

all: C

library:
	( cd lib ; $(MAKE) )

clean:
	( cd lib ; $(MAKE) clean )

purge:
	( cd lib ; $(MAKE) purge )

distclean: purge

# install ssparse
install:
	$(CP) lib/libssparse.a $(INSTALL_LIB)/libssparse.$(VERSION).a
	( cd $(INSTALL_LIB) ; ln -sf libssparse.$(VERSION).a libssparse.a )
	$(CP) include/ssparse.h $(INSTALL_INCLUDE)
	chmod 644 $(INSTALL_LIB)/libssparse*.a
	chmod 644 $(INSTALL_INCLUDE)/ssparse.h

# uninstall ssparse
uninstall:
	$(RM) $(INSTALL_LIB)/libssparse*.a
	$(RM) $(INSTALL_INCLUDE)/ssparse.h

