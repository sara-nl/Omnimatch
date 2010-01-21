#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw.h>
#include <rfftw.h>

/* routine sorts floatdata for fftw IN PLACE */
void sort4fftw(fftw_real *ext_array, float *floatdata, int Nx, int Ny, int Nz)
{
  int Nx_red, iz, ix, iy, lauf;

  Nx_red = (int) Nx/2 +1;
  lauf = 0;
  for (iz = 0; iz < Nz; iz++)
    {
      for (iy = 0; iy < Ny; iy++)
	{
	  for (ix = 0; ix < Nx; ix++)
	    {
	      ext_array[ix + 2*Nx_red * (iy + Ny * iz)] = floatdata[lauf];
	      lauf++;
	    }
	}
    }
}

/* reverse sorting after FFTW in place */
void sortback4fftw(fftw_real *ext_array, float *floatdata, int Nx, int Ny, int Nz)
{
  int Nx_red, iz, ix, iy, lauf;

  Nx_red = (int) Nx/2 +1;
  lauf = 0;
  for (iz = 0; iz < Nz; iz++)
    {
      for (iy = 0; iy < Ny; iy++)
	{
	  for (ix = 0; ix < Nx; ix++)
	    {
	      floatdata[lauf] = (float) ext_array[ix + 2*Nx_red * (iy + Ny * iz)];
	      lauf++;
	    }
	}
    }
}











