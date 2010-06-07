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

/* routine sorts floatdata for fftw IN PLACE */
void sort4fftw( fftw_real *ext_array, float *floatdata, int Nx, int Ny, int Nz ) {
  int Nx_red, iz, ix, iy, lauf;

  Nx_red = ( int ) Nx / 2 + 1;
  lauf = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        ext_array[ix + 2*Nx_red * ( iy + Ny * iz )] = floatdata[lauf];
        lauf++;
      }
    }
  }
}

/* reverse sorting after FFTW in place */
void sortback4fftw( fftw_real *ext_array, float *floatdata, int Nx, int Ny, int Nz ) {
  int Nx_red, iz, ix, iy, lauf;

  Nx_red = ( int ) Nx / 2 + 1;
  lauf = 0;
  for ( iz = 0; iz < Nz; iz++ ) {
    for ( iy = 0; iy < Ny; iy++ ) {
      for ( ix = 0; ix < Nx; ix++ ) {
        floatdata[lauf] = ( float ) ext_array[ix + 2*Nx_red * ( iy + Ny * iz )];
        lauf++;
      }
    }
  }
}











