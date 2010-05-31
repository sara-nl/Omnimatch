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

void cross( float *volinout, int Vx_max )
/* procedure exchanges quadrants
only for cubic volumes of even dimensions
FF 01/13/02
*/
{
  int lauf, i, j, k, tmpx, tmpy, tmpz;
  float tempval;

  lauf = 0;
  for ( i = 0; i < Vx_max / 2; i++ ) {
    for ( j = 0 ; j < Vx_max; j++ ) {
      for ( k = 0; k < Vx_max; k++ ) {
        if ( i < ( Vx_max / 2 ) ) tmpx = i + Vx_max / 2;
        else tmpx = i - Vx_max / 2;
        if ( j < Vx_max / 2 ) tmpy = j + Vx_max / 2;
        else tmpy = j - Vx_max / 2;
        if ( k < Vx_max / 2 ) tmpz = k + Vx_max / 2;
        else tmpz = k - Vx_max / 2;
        //printf ("tmpz: %d ... \n", tmpz);fflush(stdout);
        tempval = volinout[tmpx*Vx_max*Vx_max+tmpy*Vx_max+tmpz];
        volinout[tmpx*Vx_max*Vx_max+tmpy*Vx_max+tmpz] = volinout[lauf];
        volinout[lauf] = tempval;
        lauf++;
      }
    }
  }
}
