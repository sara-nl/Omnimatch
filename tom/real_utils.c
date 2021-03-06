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

#include <tom.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* routine symmetrizes volume along z-axis  */
void symref( float *volume, int nfold, int Nx, int Ny, int Nz ) {
  float dphi, phi, psi, the, zero __attribute__( ( unused ) );
  float *tmpvol, *tmpvol2, *symvol;
  int iz, iy, ix, irot, irun;

  printf( "  nfold = %i Nx = %i Ny = %i Nz = %i \n", nfold , Nx, Ny, Nz );
  /* check nfold */
  if ( nfold < 2 ) {
    printf( " Stop in symref - nfold = %i \n", nfold );fflush ( stdout );
    exit ( 48 );
  }
  /* allocate memory for tmpvol */
  if ( ! ( tmpvol = ( float * ) malloc ( sizeof ( float ) * Nx * Ny * Nz ) ) ) {
    printf ( " Exit in symref ... \n" );
    printf ( " Memory allocation  failure in tmpvol !!!\n" );fflush ( stdout );
    exit ( 49 );
  }
  tmpvol2 = ( float * ) malloc ( sizeof ( float ) * Nx * Ny * Nz );
  if ( ! ( symvol = ( float * ) malloc ( sizeof ( float ) * Nx * Ny * Nz ) ) ) {
    printf ( " Exit in symref ... \n" );
    printf ( " Memory allocation  failure in symvol !!!\n" );fflush ( stdout );
    exit ( 50 );
  }
  irun = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        symvol[irun] = volume[irun];
        tmpvol2[irun] = volume[irun];
        irun++;
      }
    }
  }
  psi = 0.0;the = 0.0;
  dphi = 360.0 / ( ( float ) nfold );
  for ( irot = 2; irot < nfold + 1; irot++ ) {
    phi = dphi * ( ( float ) ( irot - 1 ) );
    //phi = 18.0;
    printf( "  nfold = %i Nx = %i Ny = %i Nz = %i \n", nfold , Nx, Ny, Nz );
    tom_rotate3d( tmpvol, tmpvol2, phi, psi, the, Nx, Ny, Nz );
    printf( "rotated by phi = %f \n", phi );
    // add rotated vol
    irun = 0;
    for ( iz = 0; iz < Nz; iz++ ) {
      for ( iy = 0; iy < Ny; iy++ ) {
        for ( ix = 0; ix < Nx; ix++ ) {
          symvol[irun] = symvol[irun] + tmpvol[irun];
          irun++;
        }
      }
    }
  }
  /* divide by nfold  */
  irun = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        //volume[irun] = tmpvol[irun]/((float) nfold);
        volume[irun] = symvol[irun] / ( ( float ) nfold );
        irun++;
      }
    }
  }
}


/* routine calculates mean and variance of a volume */
float variance( float *volume, int Nx, int Ny, int Nz ) {
  int ix, iy, iz, irun;
  float mean, sq, fdims, var;

  fdims = ( ( float ) Nx * Ny * Nz );
  mean = 0;
  sq = 0;
  irun = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        mean = volume[irun] + mean;
        sq = volume[irun] * volume[irun] + sq;
        irun++;
      }
    }
  }
  mean = mean / fdims;
  sq = sq - fdims * mean * mean;
  var = sqrt( sq );
  return var;
}


/* routine normalizes volume according to variance */
void norm( float *volume, int Nx, int Ny, int Nz ) {
  int ix, iy, iz, irun;
  float mean, sq, fdims, var;

  fdims = ( ( float ) Nx * Ny * Nz );
  mean = 0;
  sq = 0;
  irun = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        mean = volume[irun] + mean;
        sq = volume[irun] * volume[irun] + sq;
        irun++;
      }
    }
  }
  mean = mean / fdims;
  sq = sq - fdims * mean * mean;
  var = sqrt( sq );
  /*  subtract mean and divide by varinace */
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        volume[irun] = volume[irun] - mean;
        volume[irun] = volume[irun] / var;
        irun++;
      }
    }
  }
}

/* limits volume to limits LOW and HI  */
void limit( float *volume, int Nx, int Ny, int Nz, float low, float hi ) {
  int ix, iy, iz, irun;

  irun = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        if ( volume[irun] < low ) volume[irun] = low;
        if ( volume[irun] > hi ) volume[irun] = hi;
        irun++;
      }
    }
  }
}


/* limits volume to limits LOW and HI - lower values are set to ZERO */
void limitz( float *volume, int Nx, int Ny, int Nz, float low, float hi ) {
  int ix, iy, iz, irun;
  float zero;

  zero = 0.0;
  irun = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        if ( volume[irun] < low ) volume[irun] = zero;
        if ( volume[irun] > hi ) volume[irun] = hi;
        irun++;
      }
    }
  }
}


/* rotates point with coordinates (x,y,z) by Euler angles phi, psi, the */
/*void pointrotate(float x, float y, float z, float phi, float psi, float theta)
{
}*/
