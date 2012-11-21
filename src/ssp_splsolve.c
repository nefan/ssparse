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
/* solve Lx=b(:,k), where L is lower triangular with unit diagonal */
CS_INT ssp_splsolve (cs *G, const cs *B, CS_INT k, CS_INT *reach, CS_ENTRY *x, const CS_INT *pinv)
{
    CS_INT j, J, p, q, px, n, *Gp, *Gi, *Bp, *Bi ;
    CS_ENTRY *Gx, *Bx ;
    CS_INT top = reach[0];
    CS_INT *xi = reach+1;

    if (!CS_CSC (G) || !CS_CSC (B) || !xi || !x) return (-1) ;
    Gp = G->p ; Gi = G->i ; Gx = G->x ; n = G->n ;
    Bp = B->p ; Bi = B->i ; Bx = B->x ;
    // top = cs_reach (G, B, k, xi, pinv) ;        /* xi[top..n-1]=Reach(B(:,k)) */
    
    for (p = top ; p < n ; p++) x [xi [p]] = 0 ;    /* clear x */
    for (p = Bp [k] ; p < Bp [k+1] ; p++) x [Bi [p]] = Bx [p] ; /* scatter B */
    for (px = top ; px < n ; px++)
    {
        j = xi [px] ;                               /* x(j) is nonzero */
        J = pinv ? (pinv [j]) : j ;                 /* j maps to col J of G */
        if (J < 0) continue ;                       /* column J is empty */
        p = Gp [J]+1 ;            /* lo: L(j,j) 1st entry */
        q = Gp [J+1] ;
        for ( ; p < q ; p++)
        {
            x [Gi [p]] -= Gx [p] * x [j] ;          /* x(i) -= G(i,j) * x(j) */
        }
    }
    return (top) ;                                  /* return top of stack */
}
