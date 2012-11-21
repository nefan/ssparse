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
/* solve Ux=b(:,k), where U is lower triangular with possible non-unit diagonal */
CS_INT ssp_sputsolve (cs *G, const cs *B, CS_INT k, CS_INT *reach, CS_ENTRY *x, const CS_INT *pinv)
{
    CS_INT j, J, p, q, px, n, *Gp, *Gi, *Bp, *Bi ;
    CS_ENTRY *Gx, *Bx ;
    CS_INT top = reach[0];
    CS_INT *xi = reach+1;

    if (!CS_CSC (G) || !CS_CSC (B) || !xi || !x) return (-1) ;
    Gp = G->p ; Gi = G->i ; Gx = G->x ; n = G->n ;
    Bp = B->p ; Bi = B->i ; Bx = B->x ;
    
    for (p = top ; p < n ; p++) x [xi [p]] = 0 ;    /* clear x */
    for (p = Bp [k] ; p < Bp [k+1] ; p++) x [Bi [p]] = Bx [p] ; /* scatter B */
#if 0
    for (px = top ; px < n ; px++)
    {
        j = xi [px] ;                               /* x(j) is nonzero */
        J = pinv ? (pinv [j]) : j ;                 /* j maps to col J of G */
        if (J < 0) continue ;                       /* column J is empty */
        p = Gp [J] ;
        q = Gp [J+1]-1 ;
        for ( ; p < q ; p++)
        {
            //if (x [Gi [p]] != 0)
                x [j] -= CS_CONJ(Gx [p]) * x [Gi [p]] ;          /* x(j) -= G(i,j) * x(i) */
        }
        x [j] /= CS_CONJ(Gx [Gp [J+1]-1]) ;
    }
    return (top) ;                                  /* return top of stack */
#else
    CS_INT lo = 1;
    for (px = top ; px < n ; px++)
    {
        j = xi [px] ;                               /* x(j) is nonzero */
        J = pinv ? (pinv [j]) : j ;                 /* j maps to col J of G */
        if (J < 0) continue ;                       /* column J is empty */
        x [j] /= Gx [lo ? (Gp [J]) : (Gp [J+1]-1)] ;/* x(j) /= G(j,j) */
        p = lo ? (Gp [J]+1) : (Gp [J]) ;            /* lo: L(j,j) 1st entry */
        q = lo ? (Gp [J+1]) : (Gp [J+1]-1) ;        /* up: U(j,j) last entry */
        for ( ; p < q ; p++)
        {
            x [Gi [p]] -= Gx [p] * x [j] ;          /* x(i) -= G(i,j) * x(j) */
        }
    }
    return (top) ;                                  /* return top of stack */
#endif
}
