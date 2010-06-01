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

/*=================================================================
 *
 * tom_rotate3d.c       Rotates a 3D Volume corresponding to the
 *                      three Euler angles, Phi, Psi, Theta
 *
 * The calling syntax is:
 *
 *              [OUT] = tom_rotate3d(IN,PHI,PSI,THETA)
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 01.07.2002
 * By: Stephan Nickell
 * Revision: 1.00 by
 *
 *=================================================================*/

#include <tom.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define pi 3.14159265


void rotate3d( float *O, float *I, float Phi, float Psi, float The, int Ox_max, int Oy_max, int Oz_max)
{
  /* Input Volumen I, Output Volumen O, Winkel the, phi, psi, Groesse der Felder xmin...zmax */
  float rotation_matrix[4][4];
  float X[4];
  float New_pt[4];
  float MP[4];
  int New_pt_ix,New_pt_iy,New_pt_iz;
  float tmp,VX1,VX2,VY1,VY2,VZ1,VZ2;
  float tmp1,tmp2,tmp3,tmp4;
  float Phirad,Psirad,Therad;
  float Cos_phi, Sin_phi, Cos_the, Sin_the, Cos_psi, Sin_psi;
  int Lay_xy;
  int i,lauf_pt,k,j;
  int lauf;
  int laufz,laufy,laufx,laufmx,laufmy;
  
  int Flag;
  
  lauf=0;
  for (i=0;i<Oz_max;i++)
    {
      for(j=0;j<Oy_max;j++)
	{
	  for(k=0;k<Ox_max;k++)
	    {
	      O[lauf]=I[lauf];
	      lauf = lauf + 1;
	    }
	}
    }
}
