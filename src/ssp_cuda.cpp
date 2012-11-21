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

// initialize CUDA
ssp_cuda *ssp_init_cuda() {
    ssp_cuda *cudaHandle = (ssp_cuda*)malloc(sizeof(ssp_cuda));
    if (!cudaHandle) {
        fprintf(stderr,"ssp_init_cuda: cudaHandle memory allocation failed.\n");
        return NULL;
    }
    cudaHandle->cusparse_handle = 0;
    cudaHandle->cusparse_matDescr = 0;

    cusparseStatus_t status = cusparseCreate(&cudaHandle->cusparse_handle);

    if (status != CUSPARSE_STATUS_SUCCESS) {
        ssp_finalize_cuda(cudaHandle);

        fprintf(stderr,"ssp_init_cuda: cusparse initialization failed.\n");
        return NULL;
    }

    status = cusparseCreateMatDescr(&cudaHandle->cusparse_matDescr); 
    if (status != CUSPARSE_STATUS_SUCCESS) {
        ssp_finalize_cuda(cudaHandle);

        fprintf(stderr,"ssp_init_cuda: cusparse matrix setup failed.\n");
        return NULL;
    }       
    cusparseSetMatType(cudaHandle->cusparse_matDescr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(cudaHandle->cusparse_matDescr,CUSPARSE_INDEX_BASE_ZERO);


    return cudaHandle;
}

// finalize CUDA
void ssp_finalize_cuda(ssp_cuda *cudaHandle) {
    if (!cudaHandle)
        return;

    if (cudaHandle->cusparse_handle)
        cusparseDestroy(cudaHandle->cusparse_handle);
    if (cudaHandle->cusparse_matDescr) 
        cusparseDestroyMatDescr(cudaHandle->cusparse_matDescr);

    free(cudaHandle);

    cudaDeviceReset();
}
