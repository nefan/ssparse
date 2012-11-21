/*
 * This file is part of ssparse.
 *
 * Copyright (c) 2006, Timothy A. Davis.
 * File derived from:
 * CXSparse: a Concise Sparse matrix package - Extended.
 * http://www.suitesparse.com
 * 
 * Modified by Stefan Sommer (shso@elektro.dtu.dk), November 2012
 *
 * ssparse is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * ssparse is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * 
 */
 
#include "ssparse.h"

/* dot product between dense vectors x,y with x having non-zero pattern xi */
CS_ENTRY ssp_spdot (CS_INT *xi, CS_ENTRY *x, CS_ENTRY *y, CS_INT n)
{
    CS_ENTRY v = 0;
    for (int i=0; i<n; i++) {
        int j = xi[i];

        if (y[j] != 0)
            v += CS_CONJ(x[j])*y[j];
    }

    return v;
}
