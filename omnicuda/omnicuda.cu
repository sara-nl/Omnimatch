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

struct PlanList {
  struct PlanList *next;
  cufftHandle plan;
  size_t size;
};

static struct PlanList *planList = NULL;


size_t getPlanSize( cufftHandle plan ) {
  PlanList *current = planList;
  while ( current )
    if ( current->plan == plan )
      return current->size;
    else current = current->next;
  return 0;
}


bool destroyPlanSize( cufftHandle plan ) {
  PlanList **current = &planList;
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
  PlanList *record = (PlanList*)malloc( sizeof( PlanList ) );
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
  cufftDestroy(plan);
  destroyPlanSize(plan);
}


void sararfftnd_one_real_to_complex(
  sararfftnd_plan p, sarafft_real *data, sarafft_complex *out
) {
  cufftResult result = cufftExecR2C( p, data, (cufftComplex*)data );
}


void sararfftnd_one_complex_to_real(
  sararfftnd_plan p, sarafft_complex *in, sarafft_real    *out
) {

}


