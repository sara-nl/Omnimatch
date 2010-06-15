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
#include <omnicuda.h>
#include <stdlib.h> // exit


sararfftnd_plan sararfft3d_create_plan(
  int nx, int ny, int nz,
  sarafft_direction dir
) {
#ifdef USE_GPUS
  sararfftnd_plan plan;
  cufftResult result = cufftPlan3d( &plan, nx, ny, nz, dir );
  if( CUFFT_SUCCESS != result )
    exit(64); // TODO better error handling (but to do that, the caller must be rewritten)
  return plan;
#else // #ifndef USE_GPUS
  return rfftw3d_create_plan( nx, ny, nz, dir, FFTW_MEASURE | FFTW_IN_PLACE );
#endif
}


void sararfftnd_one_real_to_complex(
  sararfftnd_plan p, sarafft_real *data
) {
#ifdef USE_GPUS
  // TODO: GPU implementation
#else // #ifndef USE_GPUS
  rfftwnd_one_real_to_complex( p, data, 0 );
#endif
}


void sararfftnd_one_complex_to_real(
  sararfftnd_plan p, sarafft_complex *data
) {
#ifdef USE_GPUS
  // TODO: GPU implementation
#else // #ifndef USE_GPUS
  rfftwnd_one_complex_to_real( p, data, 0 );
#endif
}

/**
 * Destroys a previously created plan.
 * The CUDA destructor returns a result code, while the fftw2 destructor is
 * a void function. For now, the result code in the CUDA destructor is
 * ignored.
 */
void sararfftnd_destroy_plan( sararfftnd_plan plan ) {
#ifdef USE_GPUS
  cufftDestroy( plan );
#else // #ifndef USE_GPUS
  rfftwnd_destroy_plan( plan );
#endif
}

void sarafft_init() {}
