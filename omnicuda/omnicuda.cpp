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

#include "omnicuda.h"
#include <stdlib.h>
#include <cuda.h>

struct PlanList {
  struct PlanList *next;
  cufftHandle plan;
  size_t size;
};

static struct PlanList *planList = NULL;

static CUcontext cuda_context;

void sarafft_init() {
  char *OMPI_COMM_WORLD_LOCAL_RANK = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
  if (NULL == OMPI_COMM_WORLD_LOCAL_RANK) {
    printf( "OMPI_COMM_WORLD_LOCAL_RANK not set!\n" );
    exit( -1 );
  }
  int localRank = atoi(OMPI_COMM_WORLD_LOCAL_RANK);
  if ( CUDA_SUCCESS != cuInit( 0 ) ) {
    printf( "cuInit failed!\n" );
    exit( -1 );
  }
  CUdevice device;
  if ( CUDA_SUCCESS != cuDeviceGet( &device, localRank ) ) {
    printf( "cuDeviceGet failed!\n" );
    exit( -1 );
  }
  if ( CUDA_SUCCESS != cuCtxCreate( &cuda_context, CU_CTX_SCHED_YIELD, device ) ) {
    printf( "cuDeviceGet failed!\n" );
    exit( -1 );
  }
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
  cufftResult result = cufftPlan3d( &plan, nx, ny, nz, dir );
  if ( CUFFT_SUCCESS != result )
    exit( -1 ); // TODO better error handling (but to do that, the caller must be rewritten)
  setPlanSize ( plan, sizeof( sarafft_real ) * nx * ny * nz );
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
  if ( CUDA_SUCCESS != cuMemAlloc( &d_data, planSize ) ) {
    printf( "cuMemAlloc failed for plansize %li!\n", planSize );
    exit( -1 );
  }
  if ( CUDA_SUCCESS != cuMemcpyHtoD( d_data, h_data, planSize ) ) {
    printf( "cuMemcpyHtoD failed!\n" );
    exit( -1 );
  }
  if ( CUFFT_SUCCESS != cufftExecR2C( plan, ( cufftReal* )d_data, ( cufftComplex* )d_data ) ) {
    printf( "cufftExecR2C failed!\n" );
    exit( -1 );
  }
  if ( CUDA_SUCCESS != cuMemcpyDtoH( h_data, d_data, planSize ) ) {
    printf( "cuMemcpyDtoH failed!\n" );
    exit( -1 );
  }
  if ( CUDA_SUCCESS != cuMemFree( d_data ) ) {
    printf( "cuMemFree failed!\n" );
    exit( -1 );
  }
}


void sararfftnd_one_complex_to_real(
  sararfftnd_plan plan, sarafft_complex *h_data
) {
  CUdeviceptr d_data;
  size_t planSize = getPlanSize( plan );
  if ( CUDA_SUCCESS != cuMemAlloc( &d_data, planSize ) ) {
    printf( "cuMemAlloc failed for plansize %li!\n", planSize );
    exit( -1 );
  }
  if ( CUDA_SUCCESS != cuMemcpyHtoD( d_data, h_data, planSize ) ) {
    printf( "cuMemcpyHtoD failed!\n" );
    exit( -1 );
  }
  if ( CUFFT_SUCCESS != cufftExecC2R( plan, ( cufftComplex* )d_data, ( cufftReal* )d_data ) ) {
    printf( "cufftExecR2C failed!\n" );
    exit( -1 );
  }
  if ( CUDA_SUCCESS != cuMemcpyDtoH( h_data, d_data, planSize ) ) {
    printf( "cuMemcpyDtoH failed!\n" );
    exit( -1 );
  }
  if ( CUDA_SUCCESS != cuMemFree( d_data ) ) {
    printf( "cuMemFree failed!\n" );
    exit( -1 );
  }
}


