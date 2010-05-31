/**************************************************************************
 * Copyright (C) 2010 W. Baumeister, MPI BioChemistry, Martinsried, Germany
 *
 * This file is part of Omnimatch.
 *
 * Omnimatch is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Omnimatch is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Omnimatch.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

//#define  pi 3.14159265359
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw.h>
#include <rfftw.h>

/* Lowpass filter (reduced complex !) */
void lowpass( fftw_complex *vol_four, float R, float smooth, int Nx, int Ny, int Nz, fftw_real scale ) {
  int ikx, iky, ikz, Nx_red, Ny_red, Nz_red;
  float dist, scf = 0.0f;
  fftw_complex zero;

  Nx_red = ( int ) Nx / 2 + 1;
  Ny_red = ( int ) Ny / 2 + 1;
  Nz_red = ( int ) Nz / 2 + 1;
  zero.re = 0.0;zero.im = 0.0;
  smooth = smooth * smooth; /* use squared smooth for Gaussian */

  for ( ikz = 0; ikz < Nz_red; ikz++ ) {
    for ( iky = 0; iky < Ny_red; iky++ ) {
      for ( ikx = 0; ikx < Nx_red; ikx++ ) {
        /* 1st quadrant */
        dist = ( ( float ) ( ikx * ikx + iky * iky + ikz * ikz ) ) ;
        dist = sqrt( dist );
        if ( dist > R ) {
          if ( smooth > 0.0 ) {
            scf = ( dist - R ) * ( dist - R );
            scf = exp( -scf / smooth );
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re * scf * scale;
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im * scf * scale;
          } else {
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny] =  zero ;
          }
        }
        /* other quadrants */
        if ( ( dist > R ) && ( iky > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny] = zero;
          }
        }
        if ( ( dist > R ) && ( ikz > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny] = zero;
          }
        }
        if ( ( dist > R ) && ( iky > 0 ) && ( ikz > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny] = zero;
          }
        }
      }
    }
  }
}

/* Bandpass - has lower and upper restrictions   */
void bandpass( fftw_complex *vol_four, float R_down, float R_up, float smooth, int Nx, int Ny, int Nz, fftw_real scale ) {
  int ikx, iky, ikz, Nx_red, Ny_red, Nz_red;
  float dist = 0.0f, scf = 0.0f;
  fftw_complex zero;

  Nx_red = ( int ) Nx / 2 + 1;
  Ny_red = ( int ) Ny / 2 + 1;
  Nz_red = ( int ) Nz / 2 + 1;
  zero.re = 0.0;zero.im = 0.0;
  smooth = smooth * smooth; /* use squared smooth for Gaussian */

  for ( ikz = 0; ikz < Nz_red; ikz++ ) {
    for ( iky = 0; iky < Ny_red; iky++ ) {
      for ( ikx = 0; ikx < Nx_red; ikx++ ) {
        dist = ( ( float ) ( ikx * ikx + iky * iky + ikz * ikz ) ) ;
        dist = sqrt( dist );
        /*   ------ LOWPASS ------  */
        /* 1st quadrant */
        if ( dist > R_up ) {
          if ( smooth > 0.0 ) {
            scf = ( dist - R_up ) * ( dist - R_up );
            scf = exp( -scf / smooth );
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re * scf * scale;
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im * scf * scale;
          } else {
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny] =  zero ;
          }
        }
        /* other quadrants */
        if ( ( dist > R_up ) && ( iky > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny] = zero;
          }
        }
        if ( ( dist > R_up ) && ( ikz > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny] = zero;
          }
        }
        if ( ( dist > R_up ) && ( iky > 0 ) && ( ikz > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].re * ( double ) scf;
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].im * ( double ) scf;
          } else {
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny] = zero;
          }
        }
        /*  ---- HIPASS ----   */
        /* 1st quadrant */
        if ( dist < R_down ) {
          if ( smooth > 0.0 ) {
            scf = ( dist - R_down ) * ( dist - R_down );
            scf = exp( -scf / smooth );
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re * scf * scale;
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im * scf * scale;
          } else {
            vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny] =  zero ;
          }
        }
        /* other quadrants */
        if ( ( dist < R_down ) && ( iky > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*( Ny-iky ) + ikz*Nx_red*Ny] = zero;
          }
        }
        if ( ( dist < R_down ) && ( ikz > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*iky + ( Nz-ikz )*Nx_red*Ny] = zero;
          }
        }
        if ( ( dist < R_down ) && ( iky > 0 ) && ( ikz > 0 ) ) {
          if ( smooth > 0 ) {
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].re =
              vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].re * ( double ) scf * scale;
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].im =
              vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny].im * ( double ) scf * scale;
          } else {
            vol_four[ikx + Nx_red*( Ny-iky ) + ( Nz-ikz )*Nx_red*Ny] = zero;
          }
        }
      }
    }
  }
}

void correl( fftw_complex *vol_four, fftw_complex *ref_four, int Nx, int Ny, int Nz, fftw_real scale ) {
  int kx, ky, kz, Nx_red, irun, ijk;
  fftw_real Ctmp, Ctmpim;

  Nx_red = Nx / 2 + 1;
  irun = 0;
  for ( kz = 0; kz < Nz; kz++ ) {
    for ( ky = 0; ky < Ny; ky++ ) {
      for ( kx = 0; kx < Nx_red; kx++ ) {
        ijk = kx + Nx_red * ( ky + kz * Ny );
        Ctmp = ref_four[ijk].re;
        Ctmpim = ref_four[ijk].im;
        ref_four[ijk].re = ( ref_four[ijk].re * vol_four[ijk].re + ref_four[ijk].im * vol_four[ijk].im ) * scale;
        ref_four[ijk].im = ( Ctmp * vol_four[ijk].im - Ctmpim * vol_four[ijk].re ) * scale;
        irun++;
      }
    }
  }
}

void convolve( fftw_complex *vol_four, fftw_complex *ref_four, int Nx, int Ny, int Nz, fftw_real scale )
/* output of convolution is written to ref_four */
{
  int kx, ky, kz, Nx_red, irun, ijk;
  fftw_real Ctmp, Ctmpim;

  Nx_red = Nx / 2 + 1;
  irun = 0;
  for ( kz = 0; kz < Nz; kz++ ) {
    for ( ky = 0; ky < Ny; ky++ ) {
      for ( kx = 0; kx < Nx_red; kx++ ) {
        ijk = kx + Nx_red * ( ky + kz * Ny );
        Ctmp = ref_four[ijk].re;
        Ctmpim = ref_four[ijk].im;
        ref_four[ijk].re = ( ref_four[ijk].re * vol_four[ijk].re - ref_four[ijk].im * vol_four[ijk].im ) * scale;
        ref_four[ijk].im = ( Ctmp * vol_four[ijk].im + Ctmpim * vol_four[ijk].re ) * scale;
        irun++;
      }
    }
  }
}
