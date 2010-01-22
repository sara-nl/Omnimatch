//#define  pi 3.14159265359
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw.h>
#include <rfftw.h>

void shift(fftw_complex *vol_four, float dx, float dy, float dz, int Nx, int Ny, int Nz)
{
  int ikx, iky, ikz, irun, Nx_red;
  float pi;
  fftw_real phase;
  fftw_real tmp;


  pi = 3.14159265359;
  Nx_red = (int) Nx/2 +1;
  /* scale shifts according to dims */
  dx = 2*pi * dx / ((float) Nx);
  dy = 2*pi * dy / ((float) Ny);
  dz = 2*pi * dz / ((float) Nz);
  
  irun = 0;
  for (ikz = 0; ikz < Nz; ikz++)
    {
      for (iky = 0; iky < Ny; iky++)
	{
	  for (ikx = 0; ikx < Nx_red; ikx++)
	    {
	      irun = ikx + Nx_red * (iky + ikz*Ny);
	      phase = (double) (-1.0 * (ikx*dx + iky*dy + ikz*dz));
	      /* multiply phase exp(i*phase) = cos(phase) + i sin(phase) */
	      tmp = vol_four[irun].re * cos(phase) -  vol_four[irun].im * sin(phase);
	      vol_four[irun].im = vol_four[irun].im * cos(phase) +   vol_four[irun].re *  sin(phase);
	      vol_four[irun].re = tmp;
	    }
	}
    }
}