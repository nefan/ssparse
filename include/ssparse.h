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
 
#ifndef _SSPS_H
#define _SSPS_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#ifdef SSP_CUDA
#include <cuda_runtime.h>
#include "cusparse_v2.h"
#endif

#include "SuiteSparse_config.h"
#include "cs.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CS_LONG
#ifdef CS_COMPLEX
#define SSP_NAME(nm) ssp_cl ## nm
#define ssp ssp_cl
#else
#define SSP_NAME(nm) ssp_dl ## nm
#define ssp ssp_dl
#endif
#else
#ifdef CS_COMPLEX
#define SSP_NAME(nm) ssp_ci ## nm
#define ssp ssp_ci
#else
#define SSP_NAME(nm) ssp_di ## nm
#define ssp ssp_di
#endif
#endif

int ssp_di_splsolve (cs_di *L, const cs_di *B, int k, int *xi, double *x,
    const int *pinv) ;
cs_long_t ssp_dl_splsolve (cs_dl *L, const cs_dl *B, cs_long_t k, cs_long_t *xi,
    double *x, const cs_long_t *pinv) ;
int ssp_ci_splsolve (cs_ci *L, const cs_ci *B, int k, int *xi, 
    cs_complex_t *x, const int *pinv) ;
cs_long_t ssp_cl_splsolve (cs_cl *L, const cs_cl *B, cs_long_t k, cs_long_t *xi, 
    cs_complex_t *x, const cs_long_t *pinv) ;

int ssp_di_sputsolve (cs_di *L, const cs_di *B, int k, int *xi, double *x,
    const int *pinv) ;
cs_long_t ssp_dl_sputsolve (cs_dl *L, const cs_dl *B, cs_long_t k, cs_long_t *xi,
    double *x, const cs_long_t *pinv) ;
int ssp_ci_sputsolve (cs_ci *L, const cs_ci *B, int k, int *xi, 
    cs_complex_t *x, const int *pinv) ;
cs_long_t ssp_cl_sputsolve (cs_cl *L, const cs_cl *B, cs_long_t k, cs_long_t *xi, 
    cs_complex_t *x, const cs_long_t *pinv) ;

CS_ENTRY ssp_di_spdot (CS_INT *xi, CS_ENTRY *x, CS_ENTRY *y, CS_INT n) ;
CS_ENTRY ssp_dl_spdot (CS_INT *xi, CS_ENTRY *x, CS_ENTRY *y, CS_INT n) ;
CS_ENTRY ssp_ci_spdot (CS_INT *xi, CS_ENTRY *x, CS_ENTRY *y, CS_INT n) ;
CS_ENTRY ssp_cl_spdot (CS_INT *xi, CS_ENTRY *x, CS_ENTRY *y, CS_INT n) ;

#define ssp_splsolve SSP_NAME (_splsolve)
#define ssp_sputsolve SSP_NAME (_sputsolve)
#define ssp_spdot SSP_NAME (_spdot)

#ifdef SSP_CUDA
typedef struct {
    cusparseHandle_t cusparse_handle;
    cusparseMatDescr_t cusparse_matDescr;
} ssp_cuda;

ssp_cuda * ssp_init_cuda(); // initialize CUDA
void ssp_finalize_cuda(ssp_cuda *cudaHandle); // close CUDA

int ssp_spdot_cuda_init (int *xpHost, int *xiHost, CS_INT n, CS_INT m, CS_INT cols, int **xpDev, int **xiDev, cuDoubleComplex **xDev, cuDoubleComplex **yDev, cuDoubleComplex **resDev, ssp_cuda *cudaHandle);
void ssp_spdot_cuda_finalize (int **xpDev, int **xiDev, cuDoubleComplex **xDev, cuDoubleComplex **yDev, cuDoubleComplex **resDev);
int ssp_spdot_cuda (cuDoubleComplex *xHost, cuDoubleComplex *yHost, cuDoubleComplex *resHost, cuDoubleComplex alpha, cuDoubleComplex beta, int *xpDev, int *xiDev,
        cuDoubleComplex *xDev, cuDoubleComplex *yDev, cuDoubleComplex *resDev, CS_INT n, CS_INT m, CS_INT cols, CS_INT nnzx, ssp_cuda *cudaHandle);

#endif


#ifdef __cplusplus
}
#endif
#endif
