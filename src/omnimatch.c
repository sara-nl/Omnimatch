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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tom.h>
#include <mpi.h>

MPI_Status status;
MPI_Request request, request1, request2, request3, request4;

int myrank, mysize, tag = 99;

void tack()
{
  time_t lt = time( NULL );
  struct tm *stop = localtime( &lt );
  printf( "Time : %i:%i:%i\n", stop->tm_hour, stop->tm_min, stop->tm_sec );
}

int main ( int argc, char *argv[] ) {
  struct em_file inputdata1;
  struct em_file inputdata2;
  struct em_file inputdata3;
  struct em_file inputdata4;
  struct em_file outputdata;

  sarafft_real *Vol_tmpl_sort, *Volume, *e3 __attribute__ ( ( unused ) ), *PointCorr, *sqconv;
  sarafft_complex *C3, *PointVolume, *PointSq;
  sararfftnd_plan p3, pi3, r3, ri3;
  sarafft_real scale;

  struct tm *zeit __attribute__ ( ( unused ) );
  struct tm start;
  char name[200];
  int Rx_max, Ry_max, Rz_max;
  int Rx_min, Ry_min, Rz_min;
  int Vx_min, Vy_min, Vz_min;
  int Vx_max, Vy_max, Vz_max;
  float Phi, Psi, Theta, winkel_lauf, nvox;
  float *Rot_tmpl, *Vol_tmpl, *Rot_mask;
  int i, j, k, tmpx __attribute__ ( ( unused ) ), tmpy __attribute__ ( ( unused ) ), tmpz __attribute__ ( ( unused ) ), lauf_pe, ksub __attribute__ ( ( unused ) );
  int ijk __attribute__ ( ( unused ) );
  int lauf;
  float max __attribute__ ( ( unused ) ), eps;
  time_t lt __attribute__ ( ( unused ) );
  float Ctmp __attribute__ ( ( unused ) ), Ctmpim __attribute__ ( ( unused ) ), Dtmp __attribute__ ( ( unused ) ), Dtmpim __attribute__ ( ( unused ) );
  int dim_fft;
  int sub[3], range[3], range_sub[3], subc[3], offset[3], dimarray[3];
  int FullVolume_dims[3];
  int nr[3];
  int area[3];

  /* MPI Variablen */
  int winkel_max, winkel_min;
  int winkel_max_pe, winkel_min_pe;
  int winkel_step_pe;
  int Phi_max, Psi_max, Theta_max;
  int Phi_min, Psi_min, Theta_min;
  int Phi_step, Psi_step, Theta_step;
  int Theta_winkel_start __attribute__ ( ( unused ) ), Psi_winkel_start __attribute__ ( ( unused ) ), Phi_winkel_start __attribute__ ( ( unused ) );
  int Theta_winkel_nr, Psi_winkel_nr, Phi_winkel_nr;
  int Theta_winkel_end __attribute__ ( ( unused ) ), Psi_winkel_end __attribute__ ( ( unused ) ), Phi_winkel_end __attribute__ ( ( unused ) );
  int Theta_steps, Psi_steps, Phi_steps;
  float Theta_winkel_rest_nr, Psi_winkel_rest_nr, Phi_winkel_rest_nr __attribute__ ( ( unused ) );
  float tempccf;
  float *Ergebnis, *conv;
  float cycles;
  int cycle;

  /* MPI Variablen Ende*/

  MPI_Init ( &argc, &argv );
  MPI_Comm_size ( MPI_COMM_WORLD, &mysize );
  MPI_Comm_rank ( MPI_COMM_WORLD, &myrank );
  sarafft_init();

  if ( argc < 15 ) {
    printf ( "\n\n" );
    printf ( " 'OMNIMATCH' computes a matched filter of tomograms \n" );
    printf ( " and templates of arbitrary geometry.\n" );
    printf ( " In detail, a locally normalized constrained cross-correlation\n" );
    printf ( " is computed (further development of CCF defined in \n" );
    printf ( "   Frangakis et al PNAS 2002. \n" );
    printf ( " ######################################################\n" );
    printf ( " All files in EM-V-int4 format !!!\n\n" );
    printf ( "PARAMETERS\n" );
    printf ( "Input: \n" );
    printf ( " no_cpus      number of CPUs to be used for MPI processing \n" );
    printf ( " Volume       Tomogram to be searched\n" );
    printf ( " Template     template of object - smaller dimension...\n" );
    printf ( " mask         (spherical) mask for local correlation - same dim as template\n" );
    printf ( " psf          pointspread function - same dim as template\n" );
    printf ( " angles       Euler angles that determine angular range:  \n" );
    printf ( "              Phi_min Phi_max Phi_step Psi_min Psi_max Psi_step The_min The_max The_step\n" );
    printf ( " dimfft       dimension of FFT - tomogram is divided into cubes of dimfft\n" );
    printf ( "              dimfft must be <= (smaller or equal to) the smallest dimension! \n\n" );
    printf ( "Output:\n" );
    printf ( " Out.ccf      non-normalized X-Correlation Function\n" );
    printf ( " Out.ccf.norm locally normalized X-Correlation Function\n" );
    printf ( " Out.ang      File containing the corresponding angles. The index of the correlated angle\n" );
    printf ( "              is stored -> keep your angular parameters to reconstruct the angles!\n" );
    printf ( "              loop order: inner=phi, outer: theta\n\n" );
    printf ( " ######################################################\n" );
    printf ( " usage: omnimatch.exe Volume Template Out ...\n" );
    printf ( "         ... Phi_min Phi_max Phi_step Psi_min Psi_max Psi_step The_min The_max The_step\n" );
    printf ( "    ... psf mask-file dim_of_fft\n\n" );
    printf ( " the total number of angles should be modulo\n" );
    printf ( " of used processors for highest computational efficiency \n\n" );
    printf ( " Linux:  1.'lamboot' to start MPI\n" );
    printf ( "              2.'mpirun -np 2 omnimatch.exe Volume Templ Out 30 180 30 30 180 30 30 180 30 Poinspread-function mask-file 256'\n\n" );
    printf ( " ######################################################\n" );
    printf ( "   Attention: \n" );
    printf ( "    Nomenclature of Euler angles is according to the EM program\n" );
    printf ( "    In most other programs (SPIDER, EMAN) Phi and Psi are \n" );
    printf ( "    exchanged compared to the EM convention. \n" );
    printf ( " ######################################################\n" );
    printf ( " In this version asymmetric masks can be used ! \n" );
    printf ( " last revision  ,  04/01/05, Friedrich Foerster\n" );
    printf ( "   - released modulu requirement " );
    printf ( " \n\n" );
    exit ( 1 );
  }

  /* Dimensionen auslesen */
  // Dimension of fft
  dim_fft = atoi ( argv[15] );
  nr[0] = 1;
  nr[1] = 1;
  nr[2] = 1;
  area[0] = dim_fft;
  area[1] = dim_fft;
  area[2] = dim_fft;
  read_em_header ( argv[1], &inputdata1 ); /* Searchvolume */
  read_em ( argv[2], &inputdata2 ); /* Template */
  FullVolume_dims[0] = inputdata1.dims[0];
  FullVolume_dims[1] = inputdata1.dims[1];
  FullVolume_dims[2] = inputdata1.dims[2];
  if ( myrank == 0 )
    printf ( " Dims tomogram: %i , %i , %i  -  dim FFT %i \n", FullVolume_dims[0], FullVolume_dims[1], FullVolume_dims[2], dim_fft );
  if ( dim_fft > FullVolume_dims[2] || dim_fft > FullVolume_dims[1] || dim_fft > FullVolume_dims[0] ) {
    if ( myrank == 0 )
      printf ( "dimfft greater than one of the dimensions! Choose smaller and restart! \n" );
    MPI_Finalize();
    exit ( 2 );
  }

  Rx_min = 1;
  Ry_min = 1;
  Rz_min = 1;
  Rx_max = ( inputdata2.dims[0] );
  Ry_max = ( inputdata2.dims[1] );
  Rz_max = ( inputdata2.dims[2] );
  Vx_min = 1;
  Vy_min = 1;
  Vz_min = 1;
  Vx_max = dim_fft;
  Vy_max = dim_fft;
  Vz_max = dim_fft;
  
    printf( "cufftPlan3d() about to start!\n" );
    fflush ( stdout );
  p3 = sararfft3d_create_plan ( Vx_max, Vy_max, Vz_max, SARAFFT_REAL_TO_COMPLEX );
  pi3 = sararfft3d_create_plan ( Vx_max, Vy_max, Vz_max, SARAFFT_COMPLEX_TO_REAL );
  r3 = sararfft3d_create_plan ( Rx_max, Rx_max, Rx_max, SARAFFT_REAL_TO_COMPLEX );
  ri3 = sararfft3d_create_plan ( Rx_max, Rx_max, Rx_max, SARAFFT_COMPLEX_TO_REAL );
  // TODO error handling for plan creation.
  if ( myrank == 0 ) {
    printf ( "Plans for FFTW created \n" );fflush ( stdout );
  }
  Volume = ( sarafft_real * ) calloc ( Vz_max * Vy_max * 2 * ( Vx_max / 2 + 1 ), sizeof ( sarafft_real ) );
  Rot_tmpl = ( float * ) malloc ( sizeof ( float ) * Rx_max * Ry_max * Rz_max );
  Rot_mask = ( float * ) malloc ( sizeof ( float ) * Rx_max * Ry_max * Rz_max );
  Vol_tmpl = ( float * ) malloc ( sizeof ( float ) * Vx_max * Vy_max * Vz_max );
  conv = ( float * ) malloc ( sizeof ( float ) * Vx_max * Vy_max * Vz_max );
  sqconv = ( sarafft_real * ) calloc ( Vz_max * Vy_max * 2 * ( Vx_max / 2 + 1 ), sizeof ( sarafft_real ) );
  if ( !
       ( inputdata1.floatdata =
           ( float * ) malloc ( sizeof ( float ) * Vx_max * Vy_max * Vz_max ) ) ) {
    printf ( "Memory allocation  failure in inputdata1.floatdata!!!" );
    fflush ( stdout );
    MPI_Finalize();
    exit ( 3 );
  }
  if ( !
       ( outputdata.floatdata =
           ( float * ) malloc ( sizeof ( float ) * Vx_max * Vy_max * Vz_max ) ) ) {
    printf ( "Memory allocation  failure in outputdata.floatdata!!!" );
    fflush ( stdout );
    MPI_Finalize();
    exit ( 4 );
  }

  if ( !
       ( Vol_tmpl_sort = ( sarafft_real * ) calloc ( Vz_max * Vy_max * 2 * ( Vx_max / 2 + 1 ), sizeof ( sarafft_real ) ) ) ) {
    printf ( "Memory allocation  failure in Volume_tmpl_sort!!!" );
    printf ( "Nx = %i, Ny = %i, Nz = %i, bytes = %li \n", 2 * ( Vx_max / 2 + 1 ), Vy_max, Vz_max, sizeof ( sarafft_real ) );
    fflush ( stdout );
    MPI_Finalize();
    exit ( 5 );
  }
  Ergebnis = ( float * ) calloc ( Vz_max * Vy_max * Vx_max, sizeof ( float ) );
  /* Winkelraum */
  Phi_min = atof ( argv[4] );
  Phi_max = atof ( argv[5] );
  Phi_step = atof ( argv[6] );
  Psi_min = atof ( argv[7] );
  Psi_max = atof ( argv[8] );
  Psi_step = atof ( argv[9] );
  Theta_min = atof ( argv[10] );
  Theta_max = atof ( argv[11] );
  Theta_step = atof ( argv[12] );
  /* Pointspread Function*/
  read_em ( argv[13], &inputdata3 );
  /* mask function */
  read_em ( argv[14], &inputdata4 );
  Phi_steps = ( Phi_max - Phi_min ) / Phi_step + 1;
  Psi_steps = ( Psi_max - Psi_min ) / Psi_step + 1;
  Theta_steps = ( Theta_max - Theta_min ) / Theta_step + 1;
  winkel_max = Phi_steps * Psi_steps * Theta_steps;
  winkel_min = 0;
  range[0] = dim_fft - 1;
  range[1] = dim_fft - 1;
  range[2] = dim_fft - 1;
  range_sub[0] = range[0] - Rx_max;
  range_sub[1] = range[1] - Rx_max;
  range_sub[2] = range[2] - Rx_max;
  sub[0] = 1;
  sub[1] = 1;
  sub[2] = 1;
  cycles = ( int ) ( FullVolume_dims[2] / ( dim_fft - Rx_max ) + 0.5 );/* cycles in z */
  cycles = ( int ) ( FullVolume_dims[1] / ( dim_fft - Rx_max ) + 0.5 ) * cycles;
  cycles = ( int ) ( FullVolume_dims[0] / ( dim_fft - Rx_max ) + 0.5 ) * cycles;
  cycle = 0;
  if ( myrank == 0 ) {
    printf ( "\n OMNIMATCH knocks the stuffing out of you ... " );tack ( &start );fflush ( stdout );
    /* prepare Output */
    strcpy ( name, argv[3] );
    strcat ( name, ".ccf" );
    printf ( "\nCreate outputfile: %s ... \n", name );fflush ( stdout );
    create_em ( name, FullVolume_dims );
    strcpy ( name, argv[3] );
    strcat ( name, ".ang" );
    printf ( "Create outputfile: %s ... \n", name );fflush ( stdout );
    create_em ( name, FullVolume_dims );
  }
  for ( sub[2] = 1; sub[2] < FullVolume_dims[2] - Rz_max;sub[2] = sub[2] + dim_fft - Rz_max ) {
    if ( myrank == 0 ) {
      tack ( &start );
      printf ( "%f%%..", ( float ) ( cycle / cycles * 100 ) );
      fflush ( stdout );
    }

    for ( sub[1] = 1; sub[1] < FullVolume_dims[1] - Ry_max;sub[1] = sub[1] + dim_fft - Ry_max ) {
      for ( sub[0] = 1; sub[0] < FullVolume_dims[0] - Rx_max;sub[0] = sub[0] + dim_fft - Rx_max ) {
        cycle = cycle + 1;
        subc[0] = sub[0];
        subc[1] = sub[1];
        subc[2] = sub[2];
        if ( sub[2] + range[2] > FullVolume_dims[2] ) subc[2] = FullVolume_dims[2] - range[2];  /* we are at the corner ?!*/
        if ( sub[1] + range[1] > FullVolume_dims[1] ) subc[1] = FullVolume_dims[1] - range[1];  /* we are at the corner ?!*/
        if ( sub[0] + range[0] > FullVolume_dims[0] ) subc[0] = FullVolume_dims[0] - range[0];  /* we are at the corner ?!*/
        read_em_subregion ( argv[1], &inputdata1, subc, range );
        read_em_subregion ( argv[1], &outputdata, subc, range );
        /* Umsortieren der Daten */
        lauf = 0;
        for ( k = 0; k < Vz_max; k++ ) {
          for ( j = 0; j < Vy_max; j++ ) {
            for ( i = 0; i < Vx_max; i++ ) {
              /* square - needed for normalization */
              sqconv[i + 2 * ( Vx_max / 2 + 1 ) * ( j + Vy_max * k ) ] = inputdata1.floatdata[lauf] * inputdata1.floatdata[lauf];
              Volume[i + 2 * ( Vx_max / 2 + 1 ) * ( j + Vy_max * k ) ] = inputdata1.floatdata[lauf];
              inputdata1.floatdata[lauf] = -1.0; /* kleine Zahl wg Max-Op , hier kommen die CCFs rein*/
              outputdata.floatdata[lauf] = -1.0; /* hier kommen die Winkel rein*/
              lauf++;
            }
          }
        }
        sararfftnd_one_real_to_complex ( p3, &Volume[0]); /* einmalige fft von Suchvolumen */
        sararfftnd_one_real_to_complex ( p3, &sqconv[0]); /* FFT of square*/
        if ( myrank < mysize - 1 ) {
          winkel_step_pe = ( int ) winkel_max / mysize;
          winkel_min_pe = myrank * winkel_step_pe;
        } else {
          winkel_step_pe = ( int ) winkel_max / mysize +
                           winkel_max - ( ( int ) winkel_max / mysize ) * mysize; /* last processor needs to do more */
          winkel_min_pe = myrank * ( ( int ) winkel_max / mysize );
        }

        if ( myrank < mysize )
          winkel_max_pe = winkel_min_pe + winkel_step_pe;
        else
          winkel_max_pe = winkel_max;/* last processor */
        Theta_winkel_nr = ( int ) winkel_min_pe / ( Psi_steps * Phi_steps );
        Theta_winkel_rest_nr =
          winkel_min_pe - Theta_winkel_nr * ( Psi_steps * Phi_steps );
        Psi_winkel_nr = ( int ) Theta_winkel_rest_nr / ( Phi_steps );
        Psi_winkel_rest_nr = Theta_winkel_rest_nr - Psi_winkel_nr * ( Phi_steps );
        Phi_winkel_nr = ( int ) Psi_winkel_rest_nr;
        Theta = Theta_winkel_nr * Theta_step + Theta_min;
        Phi = Phi_winkel_nr * Phi_step + Phi_min - Phi_step;
        Psi = Psi_winkel_nr * Psi_step + Psi_min;
        eps = 0.001;
        nvox = ( float ) 0.0;
        //Friedrich -> Zaehlung der voxels
        /*nvox = sumvoxel(Rx_max, Ry_max, Rz_max, &inputdata4.floatdata[0]);
        printf("main2: n = %f \n",nvox);fflush(stdout);*/
        eps = 0.001;
        for ( winkel_lauf = winkel_min_pe; winkel_lauf < winkel_max_pe;winkel_lauf++ ) {
          if ( Phi < Phi_max )
            Phi = Phi + Phi_step;
          else {
            Phi = Phi_min;
            Psi = Psi + Psi_step;
          }
          if ( Psi > Psi_max ) {
            Psi = Psi_min;
            Theta = Theta + Theta_step;
          }
          tom_rotate3d ( &Rot_tmpl[0], &inputdata2.floatdata[0], Phi, Psi, Theta, Rx_max, Ry_max, Rz_max );
          tom_rotate3d ( &Rot_mask[0], &inputdata4.floatdata[0], Phi, Psi, Theta, Rx_max, Ry_max, Rz_max );
          /*calculate Ref variance - ref is normalized in subroutine */
          nvox = prepref ( Rx_min, Rx_max, eps, &Rot_tmpl[0], &inputdata3.floatdata[0], &Rot_mask[0], r3, ri3 );
          pastes ( &Rot_tmpl[0], &Vol_tmpl[0], 1, 1, 1, Rx_max, Ry_max, Rz_max, Vx_max );
          scale = 1.0 / ( ( double ) Vx_max * ( double ) Vy_max * ( double ) Vz_max );
          sort4fftw ( &Vol_tmpl_sort[0], &Vol_tmpl[0], Vx_max, Vy_max, Vz_max );
          sararfftnd_one_real_to_complex ( p3, &Vol_tmpl_sort[0] );
          PointVolume = ( sarafft_complex * ) & Volume[0];
          C3 = ( sarafft_complex * ) & Vol_tmpl_sort[0];
          /* Correlation */
          correl ( &PointVolume[0], &C3[0], Vx_max, Vy_max, Vz_max, scale );
          /* back to real space */
          sararfftnd_one_complex_to_real ( pi3, &C3[0] );
          PointCorr = ( sarafft_real * ) & C3[0];
          /* reorder data */
          sortback4fftw ( &PointCorr[0], &Ergebnis[0], Vx_max, Vy_max, Vz_max );
          // cross
          cross ( &Ergebnis[0], Vx_max );
          /* ------------------- normalization ----------------------------------------- */
          /* paste mask into zero volume*/
          pastes ( &Rot_mask[0], &Vol_tmpl[0], 1, 1, 1, Rx_max, Ry_max, Rz_max, Vx_max );
          /* 1st local mean */
          sort4fftw ( &Vol_tmpl_sort[0], &Vol_tmpl[0], Vx_max, Vy_max, Vz_max );
          sararfftnd_one_real_to_complex ( p3, &Vol_tmpl_sort[0] );
          C3 = ( sarafft_complex * ) & Vol_tmpl_sort[0];
          /* Convolution of volume and mask */
          scale = 1.0 / ( ( double ) Vx_max * ( double ) Vy_max * ( double ) Vz_max );
          /*convolve( &PointVolume[0], &C3[0], Vx_max, Vy_max, Vz_max, scale);*/
          correl ( &PointVolume[0], &C3[0], Vx_max, Vy_max, Vz_max, scale );
          sararfftnd_one_complex_to_real ( pi3, &C3[0] );
          PointCorr = ( sarafft_real * ) & C3[0];
          /* reorder data (FFTW) */
          sortback4fftw ( &PointCorr[0], &conv[0], Vx_max, Vy_max, Vz_max );
          /* 2nd : convolution of square and resorting*/
          /* paste mask into zero volume*/
          pastes ( &Rot_mask[0], &Vol_tmpl[0], 1, 1, 1, Rx_max, Ry_max, Rz_max, Vx_max );
          sort4fftw ( &Vol_tmpl_sort[0], &Vol_tmpl[0], Vx_max, Vy_max, Vz_max );
          sararfftnd_one_real_to_complex ( p3, &Vol_tmpl_sort[0] );
          C3 = ( sarafft_complex * ) & Vol_tmpl_sort[0];
          PointSq = ( sarafft_complex * ) & sqconv[0];// set pointer to FFT of square
          /*convolve( &PointSq[0], &C3[0], Vx_max, Vy_max, Vz_max, scale);*/
          correl ( &PointSq[0], &C3[0], Vx_max, Vy_max, Vz_max, scale );
          sararfftnd_one_complex_to_real ( pi3, &C3[0] );
          PointCorr = ( sarafft_real * ) & C3[0];
          lauf = 0;
          for ( k = 0; k < Vz_max; k++ ) {
            for ( j = 0; j < Vy_max; j++ ) {
              for ( i = 0; i < Vx_max; i++ ) {
                /*local variance*/
                conv[lauf] = sqrt ( PointCorr[i + 2 * ( Vx_max / 2 + 1 ) * ( j + Vy_max * k ) ] - conv[lauf] * conv[lauf] / ( nvox ) ) ;
                lauf++;
              }
            }
          }
          cross ( &conv[0], Vx_max );
          /* 3rd: divide */
          lauf = 0;
          for ( k = 0 ; k < Vz_max  ; k++ ) {
            for ( j = 0; j < Vy_max; j++ ) {
              for ( i = 0; i < Vx_max; i++ ) {
                if ( conv[lauf] < eps ) {
                  tempccf = 0;
                } else {
                  tempccf = Ergebnis[lauf] / conv[lauf];/*divide CCF by variance*/
                }
                if ( inputdata1.floatdata[lauf] < tempccf ) {
                  inputdata1.floatdata[lauf] = tempccf;
                  outputdata.floatdata[lauf] = ( int ) winkel_lauf;
                }
                lauf++;
              }
            }
          }
        }                               /* Ende winkel_lauf */
        //FF
        MPI_Barrier ( MPI_COMM_WORLD );
        /* Ergebnisse einsammeln (myrank 0)*/
        if ( myrank == 0 ) {
          for ( lauf_pe = 1; lauf_pe < mysize; lauf_pe++ ) {
            MPI_Recv ( &Ergebnis[0], Vx_max * Vy_max * Vz_max, MPI_FLOAT, lauf_pe,
                       99, MPI_COMM_WORLD, &status );
            MPI_Recv ( &conv[0], Vx_max * Vy_max * Vz_max, MPI_FLOAT,
                       lauf_pe, 98, MPI_COMM_WORLD, &status );
            /* use conv as temporary memory for angles  */
            for ( lauf = 0; lauf < Vx_max * Vy_max * Vz_max; lauf++ ) {
              if ( inputdata1.floatdata[lauf] < Ergebnis[lauf] ) {
                inputdata1.floatdata[lauf] = Ergebnis[lauf];
                outputdata.floatdata[lauf] = conv[lauf];
              }
            }
          }
          /*Ergebnisse eingesammelt */
        }
        /* myrank > 0: Ergebnisse senden */
        else {
          MPI_Send ( inputdata1.floatdata, Vx_max * Vy_max * Vz_max, MPI_FLOAT, 0,
                     99, MPI_COMM_WORLD );
          MPI_Send ( outputdata.floatdata, Vx_max * Vy_max * Vz_max, MPI_FLOAT, 0,
                     98, MPI_COMM_WORLD );
        }
        MPI_Barrier ( MPI_COMM_WORLD );
        /* write non-normalized volume and angle idices */
        subc[0] = subc[0] + Rx_max / 2;
        subc[1] = subc[1] + Rx_max / 2;
        subc[2] = subc[2] + Rx_max / 2;
        if ( myrank == 0 ) {
          offset[0] = Rx_max / 2;
          offset[1] = Rx_max / 2;
          offset[2] = Rx_max / 2;
          dimarray[0] = dim_fft;
          dimarray[1] = dim_fft;
          dimarray[2] = dim_fft;
          strcpy ( name, argv[3] );
          strcat ( name, ".ccf" );
          write_em_subsubregion ( name, &inputdata1, subc, range_sub, offset, dimarray );
          strcpy ( name, argv[3] );
          strcat ( name, ".ang" );
          write_em_subsubregion ( name, &outputdata, subc, range_sub, offset, dimarray );
        }
        MPI_Barrier ( MPI_COMM_WORLD );
      }
    } /* these are the new brackets from the subregion_read , SN */
  }
  free ( Ergebnis );
  free ( inputdata1.floatdata );
  free ( inputdata2.floatdata );
  free ( inputdata3.floatdata );
  free ( inputdata4.floatdata );
  sararfftnd_destroy_plan ( p3 );
  sararfftnd_destroy_plan ( pi3 );
  sararfftnd_destroy_plan ( r3 );
  sararfftnd_destroy_plan ( ri3 );
  free ( Volume );
  free ( sqconv );
  free ( conv );
  free ( Rot_tmpl );
  free ( Rot_mask );
  free ( Vol_tmpl_sort );
  free ( outputdata.floatdata );
  if ( myrank == 0 ) {
    printf ( "OMNIMATCH finished succesfully. Good luck." );
    tack ( &start ); fflush ( stdout );
  }
  return MPI_Finalize();

  /* end main */
}

/* Variablen:

p3 und pi3:      Plan fuer die FFT bzw. IFFT (Volume)
r3 und ri3:      Plan fuer die FFT bzw. IFFT (Template Vol)
Volume:          Suchvolumen (sortiert fuer FFT)
PointVolume:     komplexer Pointer auf Suchvolumen fft(Volume)
sqconv:          Quadrat des Volumens (sortiert fuer FFT)
PointSq:         komplexer Pointer auf fft(sqconv)
conv:            Volumen fuer convolution Vol - mask / Temp fuer Winkel
Rot_tmpl:        Rotiertes Tmpl fuer Suche
Vol_tmpl:        Volumen mit rot. Template
Vol_tmpl_sort:   Volumen mit rot. Template sortiert fuer FFT
Vx_max...:       Dimensionen des Suchvolumens
Rx_max...:       Dimensionen des Templatevolumens
PointCorr:       realer Pointer auf Correlationsergebnis
C3:              komplexer Pointer auf fft(Vol_tmpl_sort)
Ergebnis:        max CCF, spaeter local mean
Ergebnis_winkel: Winkel Index
i,j,k,lauf:      Laufvariablen
*/














