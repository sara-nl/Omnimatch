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

/* ------------------------------------countvoxel----------------------------------
routine counts the non-zero voxels of a volume (cubic!) (e.g. a mask)
FF 05/02
*/
int countvoxel ( int dim, float *input, float eps ) {
  int lauf, i, j, k;
  int n;
  float temp;

  eps = 0.000001;
  //  printf("   Dim of mask %d \n", dim);fflush(stdout);
  //  printf("   eps (in countvoxel) %f \n", eps);fflush(stdout);

  n = 0;
  lauf = 0;
  for ( i = 0; i < dim; ++i ) {
    for ( j = 0; j < dim; ++j ) {
      for ( k = 0; k < dim; ++k ) {
        //     printf("lauf %d \n", n);
        temp = input[lauf];
        // printf("temp %d \n", temp);
        if ( temp > eps ) {
          n++;
        }
        lauf++;
      }
    }
  }
  //    printf(" N (out) %d \n", n);
  return( n );
}
/* ------------------------------------sumvoxel----------------------------------
routine sums the non-zero voxels of a volume (cubic!) (e.g. a mask)
FF 11/19/03
*/

float sumvoxel ( int dimx, int dimy, int dimz, float *input ) {
  float fn;
  int lauf;
  float eps;

  eps = ( float ) 0.0000001;
  fn = ( float ) 0.0;
  for ( lauf = 0; lauf < dimx*dimy*dimz; lauf++ ) {
    if ( input[lauf] > eps ) {
      fn = fn + input[lauf];
    }
  }
  return fn;
}



/* ------------------------------------energizer-----------------------------------
routine calculates rms and of a rotated reference for a given (constant ) pointspread function
FF 05/02
in this version the wedge (fft of inputdata) is multiplied on an input-volume Rot_tmpl,
furthermore a mask (same dimension as Rot_tmpl!) is multiplied and the mean is subtracted
THUS Rot_tmpl is overwritten!
FF 07/02
*/
float energizer
( int Rx_min, // unused?
  int Rx_max,
  float n,
  float *Rot_tmpl,  // reference to be manipulated
  float *inputdata, // pointspread-function
  float *mask,      // mask
  sararfftnd_plan r3,
  sararfftnd_plan ri3 ) {
  Rx_min += 0; // Suppress warning for unused parameter
  float rms_wedge;
  sarafft_real  *wedge, *Rot_tmpl_ft, *Rreal;
  sarafft_complex *W3, *R3;
  int i, j, k;
  int ijk;
  int lauf;
  float Rtmp, Rtmpim, Wtmp, Wtmpim, scale, SUM, SQSUM;

  //printf("Rx_min: %d, Rx_max: %d, n: %d\n", Rx_min, Rx_max, n);
  //  printf(" N (in energizer) %d \n", n);
  Rot_tmpl_ft = ( sarafft_real * ) malloc ( sizeof ( sarafft_real ) * Rx_max * Rx_max * 2 *
                                         ( Rx_max / 2 + 1 ) );
  wedge = ( sarafft_real * ) malloc ( sizeof ( sarafft_real ) * Rx_max * Rx_max * 2 *
                                   ( Rx_max / 2 + 1 ) );
  lauf = 0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; ++j ) {
      for ( k = 0; k < Rx_max; ++k ) {
        wedge[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] = inputdata[lauf];
        //     wedge[lauf] = inputdata2.floatdata[lauf];
        lauf++;
      }
    }
  }
  lauf = 0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; ++j ) {
      for ( k = 0; k < Rx_max; ++k ) {
        Rot_tmpl_ft[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] = Rot_tmpl[lauf];
        lauf++;
      }
    }
  }
  sararfftnd_one_real_to_complex ( r3, &Rot_tmpl_ft[0] );

  sararfftnd_one_real_to_complex ( r3, &wedge[0] );
  W3 = ( sarafft_complex * ) & wedge[0];
  R3 = ( sarafft_complex * ) & Rot_tmpl_ft[0];
  scale = 1.0 / ( Rx_max * Rx_max * Rx_max );
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; ++j ) {
      for ( k = 0; k < Rx_max / 2 + 1; ++k ) {
        ijk = k + ( Rx_max / 2 + 1 ) * ( j + i * Rx_max );
        Wtmp = c_re(W3[ijk]);
        Wtmpim = c_im(W3[ijk]);
        Rtmp = c_re(R3[ijk]);
        Rtmpim = c_im(R3[ijk]);
        c_re(R3[ijk]) = ( c_re(W3[ijk]) * c_re(R3[ijk])
                       - c_im(W3[ijk]) * c_im(R3[ijk]) ) * scale;
        c_im(R3[ijk]) = ( Wtmp * Rtmpim + Wtmpim * Rtmp ) * scale;
        /* Convolution of Image and Filter */

      }
    }
  }
  sararfftnd_one_complex_to_real ( ri3, &R3[0] );

  Rreal = ( sarafft_real * ) & R3[0];
  lauf = 0;
  SQSUM = 0;
  SUM = 0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; j++ ) {
      for ( k = 0; k < Rx_max; k++ ) {
        SQSUM = SQSUM + Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] * Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] * mask[lauf];
        SUM = SUM + Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] * mask[lauf];
        lauf++;
      }
    }
  }
  // printf("SQSUM : %f \n",SQSUM);
  //printf("SUM : %f \n",SUM);
  //printf("SUM : %d \n",n);
  SUM =  SUM / n ;
  SQSUM = SQSUM - n * SUM * SUM;
  rms_wedge = sqrt( SQSUM );
  // subtract mean within mask
  lauf = 0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; j++ ) {
      for ( k = 0; k < Rx_max; k++ ) {
        Rot_tmpl[lauf] = mask[lauf] * Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] -  mask[lauf] * SUM;
        lauf++;
      }
    }
  }
  //printf("Mean : %f \n",SUM);
  //printf("RMS : %f \n",rms_wedge);
  //printf("N: %d \n",lauf);
  free( Rot_tmpl_ft );
  free( wedge );
  return rms_wedge;
}
/* ------------------  ----------------------------------*/
float prepref( int Rx_min, // unused?
               int Rx_max,
               float eps,
               float *Rot_tmpl,  // reference to be manipulated
               float *inputdata, // pointspread-function
               float *mask,      // mask
               sararfftnd_plan r3,
               sararfftnd_plan ri3 ) {
  Rx_min += 0;
  float rms_wedge, nvox;
  sarafft_real  *wedge, *Rot_tmpl_ft, *Rreal;
  sarafft_complex *W3, *R3;
  int i, j, k;
  int ijk;
  int lauf;
  float Rtmp, Rtmpim, Wtmp, Wtmpim, scale, SUM, SQSUM;

  Rot_tmpl_ft = ( sarafft_real * ) malloc ( sizeof ( sarafft_real ) * Rx_max * Rx_max * 2 * ( Rx_max / 2 + 1 ) );
  wedge = ( sarafft_real * ) malloc ( sizeof ( sarafft_real ) * Rx_max * Rx_max * 2 * ( Rx_max / 2 + 1 ) );
  lauf = 0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; ++j ) {
      for ( k = 0; k < Rx_max; ++k ) {
        wedge[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] = inputdata[lauf];
        Rot_tmpl_ft[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] = Rot_tmpl[lauf];
        lauf++;
      }
    }
  }
  sararfftnd_one_real_to_complex ( r3, &Rot_tmpl_ft[0] );
  sararfftnd_one_real_to_complex ( r3, &wedge[0] );
  W3 = ( sarafft_complex * ) & wedge[0];
  R3 = ( sarafft_complex * ) & Rot_tmpl_ft[0];
  scale = 1.0 / ( Rx_max * Rx_max * Rx_max );
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; ++j ) {
      for ( k = 0; k < Rx_max / 2 + 1; ++k ) {
        ijk = k + ( Rx_max / 2 + 1 ) * ( j + i * Rx_max );
        Wtmp = c_re(W3[ijk]);
        Wtmpim = c_im(W3[ijk]);
        Rtmp = c_re(R3[ijk]);
        Rtmpim = c_im(R3[ijk]);
        c_re(R3[ijk]) = ( c_re(W3[ijk]) * c_re(R3[ijk])
                       - c_im(W3[ijk]) * c_im(R3[ijk]) ) * scale;
        c_im(R3[ijk]) = ( Wtmp * Rtmpim + Wtmpim * Rtmp ) * scale;
        /* Convolution (!) of Image and Filter */
      }
    }
  }
  sararfftnd_one_complex_to_real ( ri3, &R3[0] );
  Rreal = ( sarafft_real * ) & R3[0];
  lauf = 0;
  SQSUM = 0;
  SUM = 0;
  nvox = 0.0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; j++ ) {
      for ( k = 0; k < Rx_max; k++ ) {
        if ( mask[lauf] > eps ) {
          SQSUM = SQSUM + Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] * Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] * mask[lauf];
          SUM = SUM + Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )] * mask[lauf];
          nvox = mask[lauf] + nvox;
        }
        lauf++;
      }
    }
  }
  SUM =  SUM / nvox ;
  SQSUM = SQSUM - nvox * SUM * SUM;
  rms_wedge = sqrt( SQSUM );
  // subtract mean within mask and resort
  lauf = 0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; j++ ) {
      for ( k = 0; k < Rx_max; k++ ) {
        if ( mask[lauf] > eps ) {
          Rot_tmpl[lauf] = ( mask[lauf] * Rreal[k + 2 * ( Rx_max / 2 + 1 ) * ( j + Rx_max * i )]
                             -  mask[lauf] * SUM ) / rms_wedge;
        } else {
          Rot_tmpl[lauf] = 0.0;
        }
        lauf++;
      }
    }
  }
  free( Rot_tmpl_ft );
  free( wedge );
  return nvox;
}




/* -------------------------------energizer_norot------------------------------*/
/* energizer without pointspread
 FF 05/22/02          */
float energizer_norot
( int Rx_min, //unused?
  int Rx_max,
  int n,
  float *Tmpl ) {
  Rx_min += 0;
  float rms_wedge;
  int i, j, k;
  int lauf;
  float SUM, SQSUM;
  //printf("Rx_min: %d, Rx_max: %d, n: %d\n", Rx_min, Rx_max, n);
  //  printf(" N (in energizer) %d \n", n);
  lauf = 0;
  SQSUM = 0;
  SUM = 0;
  for ( i = 0; i < Rx_max; ++i ) {
    for ( j = 0; j < Rx_max; j++ ) {
      for ( k = 0; k < Rx_max; k++ ) {
        SQSUM = SQSUM + Tmpl[lauf] * Tmpl[lauf];
        SUM = SUM + Tmpl[lauf];
//           uncomment for subtraction of mean
        lauf++;
      }
    }
  }
  // printf("SQSUM : %f \n",SQSUM);
  //printf("SUM : %f \n",SUM);
  //printf("SUM : %d \n",n);
  SUM =  SUM / n ;
  SQSUM = SQSUM - n * SUM * SUM;
  rms_wedge = sqrt( SQSUM );
  return rms_wedge;
}

/*         energizer_fft
     routine calculates the local variance using FFT
     see also Roseman, A.M., Ultramicroscopy 94, p. 225-236
*/
/*float (int Rx_min,
        int Rx_max,
        int n,
        float *Rot_tmpl,  // reference to be manipulated
        float *inputdata, // pointspread-function
        float *mask,      // mask
        sararfftnd_plan r3,
        sararfftnd_plan ri3)
{
  float rms_wedge;
  sarafft_real  *wedge, *Rot_tmpl_ft, *Rreal;
  sarafft_complex *W3, *R3;
  int i, j, k;
  int ijk;
  int lauf;
  float Rtmp, Rtmpim, Wtmp, Wtmpim, scale, SUM, SQSUM;
}*/


