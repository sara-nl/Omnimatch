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

#ifndef SARAFFT_H
#define SARAFFT_H
#ifdef __cplusplus
extern "C" {
#endif



#ifdef USE_GPUS


#include <cufft.h>
// TODO: GPU declarations
typedef cufftReal sarafft_real;
typedef cufftComplex sarafft_complex;
typedef cufftHandle sararfftnd_plan;
typedef cufftType sarafft_direction;
#define SARAFFT_REAL_TO_COMPLEX CUFFT_R2C
#define SARAFFT_COMPLEX_TO_REAL CUFFT_C2R
#define c_re(f) ((f).x)
#define c_im(f) ((f).y)


#else // #ifdef USE_GPUS


#include <sfftw.h>
#include <srfftw.h>

typedef fftw_real sarafft_real;
typedef fftw_complex sarafft_complex;
typedef rfftwnd_plan sararfftnd_plan;
typedef fftw_direction sarafft_direction;
#define SARAFFT_REAL_TO_COMPLEX FFTW_FORWARD
#define SARAFFT_COMPLEX_TO_REAL FFTW_BACKWARD


#endif // #ifndef USE_GPUS


void sarafft_init();
sararfftnd_plan sararfft3d_create_plan( int nx, int ny, int nz, sarafft_direction dir );
void sararfftnd_one_real_to_complex( sararfftnd_plan p, sarafft_real    *data );
void sararfftnd_one_complex_to_real( sararfftnd_plan p, sarafft_complex *data );
void sararfftnd_destroy_plan( sararfftnd_plan plan );



#ifdef __cplusplus
} // extern "C"
#endif
#endif // #ifndef SARAFFT_H
