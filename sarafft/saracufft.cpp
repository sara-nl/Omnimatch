/**************************************************************************
 * Copyright (C) 2010 Pieter van Beek
 *
 * This file is part of SARAFFT.
 *
 * SARAFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SARAFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SARAFFT.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

#include "sarafft.h"
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

struct PlanList {
  struct PlanList *next;
  cufftHandle plan;
  size_t size;
};

static struct PlanList *planList = NULL;

static CUcontext cuda_context;

extern "C" void sarafft_init() {
  printf( "Cuda is about to be initialized!\n" );
  fflush ( stdout );
  char *OMPI_COMM_WORLD_LOCAL_RANK = getenv( "OMPI_COMM_WORLD_LOCAL_RANK" );
  if ( NULL == OMPI_COMM_WORLD_LOCAL_RANK ) {
    printf( "OMPI_COMM_WORLD_LOCAL_RANK not set!\n" );
    fflush ( stdout );
    exit( 80 );
  }
  int localRank = atoi( OMPI_COMM_WORLD_LOCAL_RANK );
  printf( "Local rank is %d\n", localRank );
  fflush ( stdout );
  if ( CUDA_SUCCESS != cuInit( 0 ) ) {
    printf( "cuInit failed!\n" );
    fflush ( stdout );
    exit( 81 );
  }
  CUdevice device;
  if ( CUDA_SUCCESS != cuDeviceGet( &device, localRank ) ) {
    printf( "cuDeviceGet failed!\n" );
    fflush ( stdout );
    exit( 82 );
  }
  if ( CUDA_SUCCESS != cuCtxCreate( &cuda_context, CU_CTX_SCHED_YIELD, device ) ) {
    printf( "cuCtxCreate failed!\n" );
    fflush ( stdout );
    exit( 83 );
  }
  printf( "Cuda was initialized successfully!\n" );
  fflush ( stdout );
}


size_t getPlanSize( cufftHandle plan ) {
  struct PlanList *current = planList;
  while ( current )
    if ( current->plan == plan )
      return current->size;
    else
      current = current->next;
  return 0;
}


bool destroyPlanSize( cufftHandle plan ) {
  struct PlanList **current = &planList;
  while ( *current )
    if ( ( *current )->plan == plan ) {
      PlanList *next = ( *current )->next;
      free( ( void* )( *current ) );
      *current = next;
      return true;
    } else current = &( ( *current )->next );
  return false;
}


void setPlanSize( cufftHandle plan, size_t size ) {
  destroyPlanSize( plan );
  struct PlanList *record = ( PlanList* )malloc( sizeof( PlanList ) );
  record->plan = plan;
  record->size = size;
  record->next = planList;
  planList = record;
}


sararfftnd_plan sararfft3d_create_plan(
  int nx, int ny, int nz, sarafft_direction dir
) {
  sararfftnd_plan plan;
  printf( "cufftPlan3d() about to start!\n" );
  fflush ( stdout );
  cufftResult result = cufftPlan3d( &plan, nx, ny, nz, dir );
  if ( CUFFT_SUCCESS != result ) {
    printf( "cufftPlan3d() failed with code %d for dir=%d\n", result, dir );
    fflush ( stdout );
    exit( 84 ); // TODO better error handling (but to do that, the caller must be rewritten)
  }
  printf( "cufftPlan3d() succeeded!\n" );
  fflush ( stdout );
  setPlanSize ( plan, sizeof( sarafft_real ) * nx * ny * ( nz + 2 ) );
  return plan;
}


void sararfftnd_destroy_plan(
  sararfftnd_plan plan
) {
  cufftDestroy( plan );
  destroyPlanSize( plan );
}


void sararfftnd_one_real_to_complex(
  sararfftnd_plan plan, sarafft_real *h_data
) {
  CUdeviceptr d_data;
  size_t planSize = getPlanSize( plan );
//   printf( "planSize = %li!\n", planSize );
//   fflush ( stdout );
  cufftResult fftResult;
  CUresult cudaResult;
  if ( CUDA_SUCCESS != cuMemAlloc( &d_data, planSize ) ) {
    printf( "cuMemAlloc failed for plansize %li!\n", planSize );
    fflush ( stdout );
    exit( 85 );
  }
  if ( CUDA_SUCCESS != cuMemcpyHtoD( d_data, h_data, planSize ) ) {
    printf( "cuMemcpyHtoD failed!\n" );
    fflush ( stdout );
    exit( 86 );
  }
//   cudaError_t cudaError = cudaGetLastError();
//   if( cudaError != cudaSuccess ) {
//     printf( "CUDA Runtime API Error reported : %s\n", cudaGetErrorString(cudaError));
//     fflush ( stdout );
//     exit( 87 );
//   } else {
//     printf( "CUDA is in good shape.\n");
//     fflush ( stdout );
//   }
  fftResult = cufftExecR2C( plan, ( cufftReal* )d_data, ( cufftComplex* )d_data );
  if ( CUFFT_SUCCESS != fftResult ) {
    printf( "cufftExecR2C failed with code %d\n", fftResult );
    fflush ( stdout );
    exit( 87 );
  }
  if ( CUDA_SUCCESS != cuMemcpyDtoH( h_data, d_data, planSize ) ) {
    printf( "cuMemcpyDtoH failed!\n" );
    fflush ( stdout );
    exit( 88 );
  }
  if ( CUDA_SUCCESS != cuMemFree( d_data ) ) {
    printf( "cuMemFree failed!\n" );
    fflush ( stdout );
    exit( 89 );
  }
}


void sararfftnd_one_complex_to_real(
  sararfftnd_plan plan, sarafft_complex *h_data
) {
  CUdeviceptr d_data;
  size_t planSize = getPlanSize( plan );
  if ( CUDA_SUCCESS != cuMemAlloc( &d_data, planSize ) ) {
    printf( "cuMemAlloc failed for plansize %li!\n", planSize );
    fflush ( stdout );
    exit( 90 );
  }
  if ( CUDA_SUCCESS != cuMemcpyHtoD( d_data, h_data, planSize ) ) {
    printf( "cuMemcpyHtoD failed!\n" );
    fflush ( stdout );
    exit( 91 );
  }
  if ( CUFFT_SUCCESS != cufftExecC2R( plan, ( cufftComplex* )d_data, ( cufftReal* )d_data ) ) {
    printf( "cufftExecR2C failed!\n" );
    fflush ( stdout );
    exit( 92 );
  }
  if ( CUDA_SUCCESS != cuMemcpyDtoH( h_data, d_data, planSize ) ) {
    printf( "cuMemcpyDtoH failed!\n" );
    fflush ( stdout );
    exit( 93 );
  }
  if ( CUDA_SUCCESS != cuMemFree( d_data ) ) {
    printf( "cuMemFree failed!\n" );
    fflush ( stdout );
    exit( 94 );
  }
}


