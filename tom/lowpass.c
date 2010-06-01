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
//#define  pi 3.14159265359
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Lowpass filter (reduced complex !) */
void lowpass(fftw_complex *vol_four, float R, float smooth, int Nx, int Ny, int Nz)
{
  int ikx, iky, ikz, Nx_red, Ny_red, Nz_red;
  float dist, scf;
  fftw_complex zero; 

  Nx_red = (int) Nx/2 +1;
  Ny_red = (int) Ny/2 +1;
  Nz_red = (int) Nz/2 +1;
  zero.re = 0.0;zero.im = 0.0;
  smooth = smooth*smooth; /* use squared smooth for Gaussian */
  
  for (ikz = 0; ikz < Nz_red; ikz++)
    {
      for (iky = 0; iky < Ny_red; iky++)
	{
	  for (ikx = 0; ikx < Nx_red; ikx++)
	    {
	      	      /* 1st quadrant */
	      dist = ((float) (ikx*ikx+iky*iky+ikz*ikz)) ;
	      dist = sqrt(dist);
	      if (dist > R) 
		{
		  if (smooth > 0.0)
		    {
		      scf = (dist-R)*(dist-R);
		      scf = exp(-scf/smooth);
		      vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re * scf;
		      vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im * scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny] =  zero ;
		    }
		}
	      /* other quadrants */
	      if ((dist > R) && (iky>0)) 
		{
		  if (smooth > 0)
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].re* (double) scf;
		      vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].im* (double) scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny] = zero;
		    }
		}
	      if ((dist > R) && (ikz>0)) 
		{
		  if (smooth > 0)
		    {
		      vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].re * (double) scf;
		      vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].im * (double) scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny] = zero;
		    }
		}
	      if ((dist > R) && (iky>0) && (ikz>0)) 
		{
		  if (smooth > 0)
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].re * (double) scf;
		      vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].im * (double) scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny] = zero;
		    }
		}
	    }
	}
    }
}

/* Bandpass - has lower and upper restrictions   */
void bandpass(fftw_complex *vol_four, float R_down, float R_up, float smooth, int Nx, int Ny, int Nz)
{
  int ikx, iky, ikz, Nx_red, Ny_red, Nz_red;
  float dist, scf;
  fftw_complex zero; 

  Nx_red = (int) Nx/2 +1;
  Ny_red = (int) Ny/2 +1;
  Nz_red = (int) Nz/2 +1;
  zero.re = 0.0;zero.im = 0.0;
  smooth = smooth*smooth; /* use squared smooth for Gaussian */
  
  for (ikz = 0; ikz < Nz_red; ikz++)
    {
      for (iky = 0; iky < Ny_red; iky++)
	{
	  for (ikx = 0; ikx < Nx_red; ikx++)
	    {
	      	      /* 1st quadrant */
	      dist = ((float) (ikx*ikx+iky*iky+ikz*ikz)) ;
	      dist = sqrt(dist);
	      if (dist > R) 
		{
		  if (smooth > 0.0)
		    {
		      scf = (dist-R)*(dist-R);
		      scf = exp(-scf/smooth);
		      vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].re * scf;
		      vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny].im * scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*iky + ikz*Nx_red*Ny] =  zero ;
		    }
		}
	      /* other quadrants */
	      if ((dist > R) && (iky>0)) 
		{
		  if (smooth > 0)
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].re* (double) scf;
		      vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny].im* (double) scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + ikz*Nx_red*Ny] = zero;
		    }
		}
	      if ((dist > R) && (ikz>0)) 
		{
		  if (smooth > 0)
		    {
		      vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].re * (double) scf;
		      vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny].im * (double) scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*iky + (Nz-ikz)*Nx_red*Ny] = zero;
		    }
		}
	      if ((dist > R) && (iky>0) && (ikz>0)) 
		{
		  if (smooth > 0)
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].re = 
			vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].re * (double) scf;
		      vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].im = 
			vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny].im * (double) scf;
		    }
		  else
		    {
		      vol_four[ikx + Nx_red*(Ny-iky) + (Nz-ikz)*Nx_red*Ny] = zero;
		    }
		}
	    }
	}
    }
}
