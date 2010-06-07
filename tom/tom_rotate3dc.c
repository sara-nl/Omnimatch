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
 * tom_rotate3d.c	Rotates a 3D Volume corresponding to the 
 *	                three Euler angles, Phi, Psi, Theta 
 *
 * The calling syntax is:
 *
 *		[OUT] = tom_rotate3d(IN,PHI,PSI,THETA)
 *
 *  Electron Tomography toolbox of the
 *  Max-Planck-Institute for Biochemistry
 *  Dept. Molecular Structural Biology
 *  82152 Martinsried, Germany
 *  http://www.biochem.mpg.de
 *
 * Last changes: 25.09.2003
 * FF
 *
 * By: Stephan Nickell
 * Revision: 1.00 by 
 *
 *=================================================================*/

#include <tom.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h>

#define pi 3.1415926535897 

/* Input Arguments */
#define	IN	prhs[0]
#define	PHI	prhs[1]
#define	PSI	prhs[2]
#define	THE	prhs[3]


/* Output Arguments */
#define	OUT	plhs[0]
  
void tom_rotate3d(
  float *O, 
  float *I, 
  float Phi, 
  float Psi, 
  float The,
  int Ox_max, 
  int Oy_max, 
  int Oz_max)
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
  
  if(Phi==0.0 && Psi==0.0 && The==0.0) {
    lauf=0;
    for (i=0;i<Ox_max;++i)
      {
	for(j=0;j<Oy_max;j++)
	  {
	    for(k=0;k<Oz_max;k++)
	      {
		O[lauf]=I[lauf];
		lauf++;
	      }
	  }
      }
  } 
  else 
    {
      Phirad=Phi*pi/180.0; /* Phi im Bogenmass */
      Psirad=Psi*pi/180.0; /* Psi im Bogenmass */
      Therad=The*pi/180.0; /* The im Bogenmass */
      
      MP[1]=((float)Ox_max/2.0);  /*Define the center of the rotation*/
      MP[2]=((float)Oy_max/2.0);
      MP[3]=((float)Oz_max/2.0);
      
      Lay_xy=(Ox_max)*(Oy_max);
      
      Cos_phi=0.0;
      Sin_phi=0.0;
      Cos_psi=0.0;
      Sin_psi=0.0;
      Cos_the=0.0;
      Sin_the=0.0;
      
      Cos_phi = cos(Phirad);
      Sin_phi = sin(Phirad);
      Cos_psi = cos(Psirad);
      Sin_psi = sin(Psirad);
      Cos_the = cos(Therad);
      Sin_the = sin(Therad);
      
      rotation_matrix[1][1]=  (Cos_psi*Cos_phi) -(Cos_the*Sin_psi*Sin_phi);
      rotation_matrix[2][1] = (Sin_psi*Cos_phi) +(Cos_the*Cos_psi*Sin_phi);
      rotation_matrix[3][1] = (Sin_the*Sin_phi);
      rotation_matrix[1][2] = (- Cos_psi*Sin_phi) -(Cos_the*Sin_psi*Cos_phi);
      rotation_matrix[2][2] = (- Sin_psi*Sin_phi) +( Cos_the*Cos_psi*Cos_phi);
      rotation_matrix[3][2] = (Sin_the*Cos_phi);
      rotation_matrix[1][3] = (Sin_the*Sin_psi);
      rotation_matrix[2][3] = (- Sin_the*Cos_psi);
      rotation_matrix[3][3] = (Cos_the); 
      
      
      lauf_pt=-1;
      lauf=0.0;
      
      
      X[3]=-MP[3]-1.0;
      for (laufz=0;laufz<Oz_max;laufz++)
	{
	  X[3]=X[3]+1.0;
	  X[2]=-MP[2]-1.0;
	  for (laufy=0;laufy<Oy_max;laufy++)
	    {
	      X[2]=X[2]+1.0;
	      X[1]=-MP[1]-1.0;
	      for (laufx=0;laufx<Ox_max;laufx++)
		{
		  lauf_pt=lauf_pt+1;
		  X[1]=X[1]+1.0;
		  Flag=0;
		  
		  
		  for (laufmx=1;laufmx<=3;laufmx++)
		    {
		      New_pt[laufmx]=MP[laufmx];
		      for(laufmy=1;laufmy<=3;laufmy++) 
			{
			  New_pt[laufmx]=New_pt[laufmx]+rotation_matrix[laufmy][laufmx]*X[laufmy];
			}
		      
		      if (New_pt[laufmx]<0.0 || New_pt[1]>(float)Ox_max-1.0 || New_pt[2]>(float)Oy_max-1.0  || New_pt[3]>(float)Oz_max-1.0)
			{
			  O[lauf_pt]=0.0;Flag=1;
		       }
		      
		    }
		  
		  
		  
		  if (Flag<1)
		    {
		      
		      tmp=New_pt[1];
		      New_pt_ix=(int) tmp;
		      VX2= tmp-(float)New_pt_ix;
		      VX1= 1.0-(float)VX2;
		      tmp=New_pt[2];
		      New_pt_iy=(int) tmp;
		      VY2=tmp-(float) New_pt_iy;
		      VY1=1.0-(float)VY2;
		      lauf=New_pt_ix+1+New_pt_iy*Ox_max;
		      
		      tmp=New_pt[3];
		      New_pt_iz=(int) tmp;
		      VZ2=tmp-(float)New_pt_iz;
		      VZ1=1.0-(float)VZ2;
		      lauf=lauf+New_pt_iz*Lay_xy-1;
		      /* bug fixed FF */
		      if ( (lauf+Lay_xy+Ox_max+1 < Ox_max*Oy_max*Oy_max) && (lauf+Ox_max+1 < Ox_max*Oy_max*Oy_max) && (lauf+Lay_xy+1 < Ox_max*Oy_max*Oy_max) && (lauf+Lay_xy+1 < Ox_max*Oy_max*Oy_max )  )
			{ 
			  tmp1=I[lauf]+(I[lauf+1]-I[lauf])*VX2;
			  tmp2=I[lauf+Ox_max]*VX1+I[lauf+Ox_max+1]*VX2;
			  tmp3=I[lauf+Lay_xy]*VX1+I[lauf+Lay_xy+1]*VX2;
			  tmp4=I[lauf+Lay_xy+Ox_max]*VX1+I[lauf+Lay_xy+Ox_max+1]*VX2;  
			  tmp1=tmp1*VY1+tmp2*VY2;
			  tmp2=tmp3*VY1+tmp4*VY2;
			  O[lauf_pt]=tmp1*VZ1+tmp2*VZ2; 
			} 
		    }
		}
	    }
	}
    }
  /* end rot */
}





