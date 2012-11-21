/*
 * This file is part of ssparse.
 *
 * Copyright (C) 2012, Technical University of Denmark
 * https://github.com/nefan/ssparse
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

/* dot product between sets of dense vectors x,y with x having non-zero pattern xi 
   CUDA accelerated version */
int ssp_spdot_cuda (cuDoubleComplex *xHost, cuDoubleComplex *yHost, cuDoubleComplex *resHost, 
        cuDoubleComplex alpha, cuDoubleComplex beta, int *xpDev, int *xiDev,
        cuDoubleComplex *xDev, cuDoubleComplex *yDev, cuDoubleComplex *resDev, 
        CS_INT n, CS_INT m, CS_INT cols, CS_INT nnzx, ssp_cuda *cudaHandle)
{
    cudaError_t stat1 = cudaMemcpy(xDev, xHost, (size_t)nnzx*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    cudaError_t stat2 = cudaMemcpy(yDev, yHost, (size_t)n*cols*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    if ((stat1 != cudaSuccess) || (stat2 != cudaSuccess)) {
        fprintf(stderr,"ssp_spdot_cuda: host->device memory copy failed.\n");
        return false;
    }

    if (beta.x != 0 && beta.y != 0) {
        stat1 = cudaMemcpy(resDev, resHost, (size_t)std::max(n,m)*cols*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
        if ((stat1 != cudaSuccess)) {
            fprintf(stderr,"ssp_spdot_cuda: host->device memory copy failed.\n");
            return false;
        }
    }

    // do inner product computation
    cusparseStatus_t status = cusparseZcsrmm(cudaHandle->cusparse_handle, 
            CUSPARSE_OPERATION_NON_TRANSPOSE, 
            m, cols, n, nnzx, &alpha, 
            cudaHandle->cusparse_matDescr, 
            xDev, xpDev, xiDev, yDev, 
            n, &beta, resDev, std::max(n,m));
    if (status != CUSPARSE_STATUS_SUCCESS) {
        fprintf(stderr,"ssp_spdot_cuda: matrix-matrix multiplication failed.\n");
        return false;
    }  
    
    stat1 = cudaMemcpy(resHost, resDev, (size_t)std::max(n,m)*cols*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    if ((stat1 != cudaSuccess)) {
        fprintf(stderr,"ssp_spdot_cuda: device->host memory copy failed. %d\n",stat1);
        return false;
    }

    return true;
}

// allocate memory etc.
int ssp_spdot_cuda_init (int *xpHost, int *xiHost, CS_INT n, CS_INT m, CS_INT cols, 
        int **xpDev, int **xiDev, 
        cuDoubleComplex **xDev, cuDoubleComplex **yDev, cuDoubleComplex **resDev, 
        ssp_cuda *cudaHandle) {

    int nnzx = xpHost[m];

    cudaError_t stat1 = cudaMalloc(xpDev, (m+1)*sizeof(int));
    cudaError_t stat2 = cudaMalloc(xiDev, nnzx*sizeof(int));
    cudaError_t stat3 = cudaMalloc(xDev, nnzx*sizeof(cuDoubleComplex));
    cudaError_t stat4 = cudaMalloc(yDev, n*cols*sizeof(cuDoubleComplex));
    cudaError_t stat5 = cudaMalloc(resDev, std::max(n,m)*cols*sizeof(cuDoubleComplex));

    if ((stat1 != cudaSuccess) || (stat2 != cudaSuccess) ||
            (stat3 != cudaSuccess) || (stat4 != cudaSuccess) ||
            (stat5 != cudaSuccess)) {

        ssp_spdot_cuda_finalize(xpDev,xiDev,xDev,yDev,resDev);

        fprintf(stderr,"ssp_spdot_cuda_init: memory allocation failed.\n");
        return false;
    }

    stat1 = cudaMemcpy(*xpDev, xpHost, (size_t)(m+1)*sizeof(int), cudaMemcpyHostToDevice);
    stat2 = cudaMemcpy(*xiDev, xiHost, (size_t)nnzx*sizeof(int), cudaMemcpyHostToDevice);
    if ((stat1 != cudaSuccess) || (stat2 != cudaSuccess)) {

        ssp_spdot_cuda_finalize(xpDev,xiDev,xDev,yDev,resDev);

        fprintf(stderr,"ssp_spdot_cuda_init: host->device memory copy failed.\n");
        return false;
    }

    return true;
}

void ssp_spdot_cuda_finalize (int **xpDev, int **xiDev, cuDoubleComplex **xDev, cuDoubleComplex **yDev, cuDoubleComplex **resDev) {
    if (*xpDev) cudaFree(*xpDev); *xpDev = 0;
    if (*xiDev) cudaFree(*xiDev); *xiDev = 0;
    if (*xDev) cudaFree(*xDev); *xDev = 0;
    if (*yDev) cudaFree(*yDev); *yDev = 0;
    if (*resDev) cudaFree(*resDev); *resDev = 0;
}
