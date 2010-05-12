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

// routine pastes volume of dimension Ox_max into bigger volume (dim Vx_max)
// volume ist pasted into the MIDDLE for EVEN and ODD dimensions of Ox-max
// odd dimension of V_x not tested
// FF 07/11/02
void pastes(float *I, float *O,int Ox_min __attribute__((unused)), int Oy_min __attribute__((unused)), int Oz_min __attribute__((unused)), int Ox_max, int Oy_max, int Oz_max,int Vx_max)
{
int laufx, laufy, laufz;
int laufp/*, laufry, laufrz*/;
int Vx_max_2,Vy_max_2 __attribute__((unused)),Vz_max_2 __attribute__((unused)),Ox_max_2,Oy_max_2,Oz_max_2;
int i,j,k;
Ox_max_2= Ox_max /2;
//Ox_max_2= Ox_max /2 + Ox_max % 2 ;
//printf("Ox: %d , Ox halbe: %d \n", Ox_max, Ox_max_2);
Oy_max_2= Oy_max /2 ;
Oz_max_2= Oz_max /2 ;
Vx_max_2=Vx_max/2-1;
laufp=0;
// routine pastes a small volume (reference) into a bigger one (O)

/* leer anlegen */
for (i=0;i<Vx_max;++i)
 	{
 	for(j=0;j<Vx_max;j++)
		{
 		for(k=0;k<Vx_max;k++)
			{
			O[laufp]=0.0;
			laufp++;
			}
		}
	}


laufp=0;
for (laufz=Vx_max_2-Ox_max_2+1;laufz<=Vx_max_2-Ox_max_2+Ox_max;laufz++)
	{
	for (laufy=Vx_max_2-Oy_max_2+1;laufy<=Vx_max_2-Oy_max_2+Oy_max;laufy++)
		{
		for (laufx=Vx_max_2-Oz_max_2+1;laufx<=Vx_max_2-Oz_max_2+Oz_max;laufx++)
			{
			O[(laufx+Vx_max*(laufy+Vx_max*laufz))]=I[laufp];
			laufp++;
			}
		}
	}


}



