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

/* Cray 'int' sind 8 bytes aber SGI 'int' sind nur 4 bytes */
/*Cray :

typedef short int int32;
*/
/*SGI :
typedef int int32;
*/
/* emfile.h schreibt und liest EM Files im int 4 stream
   format (READ V)
   analog zu EM, getestet -> exakt das gleiche Ergebnis
   letzte Aenderung 18.11.1999
   Stephan Nickell
   MPI fuer Biochemie
*/



/* already defined in tom.h:
struct em_file {

  unsigned char magic[1];
  char dummya[2];
  unsigned char type[1];

  int dims[3];

  char comment[80];
  int emdata[40];
  char dummyb[256];

  unsigned char *bytedata;
  int *intdata;
  float *floatdata;
};
*/

void read_em ( char *infile, struct em_file *inemdata ) {
  FILE *input = 0;
  long size;
//  int lauf;
  if ( ( input = fopen ( infile, "r" ) ) == 0 ) {
    printf( "could not open %s\n", infile );
    exit ( 16 );
  }


  fread ( inemdata->magic, 1, 1, input );
  fread ( inemdata->dummya, 1, 2, input );
  fread ( inemdata->type, 1, 1, input );
  fread ( inemdata->dims, 4, 3, input );
  fread ( inemdata->comment, 1, 80, input );
  fread ( inemdata->emdata, 4, 40, input );
  fread ( inemdata->dummyb, 1, 256, input );

  size = inemdata->dims[0] * inemdata->dims[1] * inemdata->dims[2];

  switch ( inemdata->type[0] ) {
  case 1: if ( ( inemdata->bytedata = ( unsigned char * )( malloc( size ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 17 ); }
    break;
  case 2: if ( ( inemdata->intdata = ( int * )( malloc( size * 2 ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 18 ); }
    if ( ( inemdata->floatdata = ( float * )( malloc( size * 4 ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 19 );}
    break;
  case 5: if ( ( inemdata->floatdata = ( float * )( malloc( size * 4 ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 20 ); }
    break;
  }


  switch ( inemdata->type[0] ) {
  case 1: fread ( inemdata->bytedata, 1, size, input );
  case 2: fread ( inemdata->intdata, 2, size, input );
  case 5: fread ( inemdata->floatdata, 4, size, input );
  }
  fclose ( input );

};

void read_em_header ( char *infile, struct em_file *inemdata ) {
  FILE *input = 0;
//     long size;
//  int lauf;
  if ( ( input = fopen ( infile, "r" ) ) == 0 ) {
    printf( "could not open %s\n in readm_em_header", infile );
    exit ( 21 );
  }

  fread ( inemdata->magic, 1, 1, input );
  fread ( inemdata->dummya, 1, 2, input );
  fread ( inemdata->type, 1, 1, input );
  fread ( inemdata->dims, 4, 3, input );
  fread ( inemdata->comment, 1, 80, input );
  fread ( inemdata->emdata, 4, 40, input );
  fread ( inemdata->dummyb, 1, 256, input );

  fclose ( input );
  /*    printf("read_em_header ... ");fflush(stdout); */
};


void read_em_subregion ( char *infile, struct em_file *inemdata, int *nr, int  *area ) {
  FILE *input = 0;
  long size;
  int lauf, ilaufx, ilaufz, laufy ;
  long int s1, s2, s3;
  int area_d[3];
  long int xy_dims, fseek_merker;
//     unsigned int laufin;
  int size_area;
  int count;

  if ( ( input = fopen ( infile, "r" ) ) == 0 ) {
    printf( "could not open %s\n", infile );
    exit ( 22 );
  }

  fread ( inemdata->magic, 1, 1, input );
  fread ( inemdata->dummya, 1, 2, input );
  fread ( inemdata->type, 1, 1, input );
  fread ( inemdata->dims, 4, 3, input );
  fread ( inemdata->comment, 1, 80, input );
  fread ( inemdata->emdata, 4, 40, input );
  fread ( inemdata->dummyb, 1, 256, input );


  size = ( nr[0] + area[0] ) * ( nr[1] + area[1] ) * ( nr[2] + area[2] );
  /*    printf("Size: %i\n",size);fflush(stdout);
  printf("Nr.: %i %i %i\n",nr[0],nr[1],nr[2]);fflush(stdout);
  printf("Area: %i %i %i\n",area[0],area[1],area[2]);fflush(stdout); */

  /*   switch(inemdata->type[0])
    {
    case 1: break;
    case 2: break;
    case 5: if ((inemdata->floatdata=(float *)(malloc(size*sizeof(float)))) == 0)
  { printf ("could not allocate memory in read_em_subregion"); exit (1); }
    break;
    }
  */
  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;


  switch ( inemdata->type[0] ) {
  case 5:
    fseek( input, 4*( nr[0] - 1 ), SEEK_CUR );
    fseek( input, 4*( inemdata->dims[0]*( nr[1] - 1 ) ), SEEK_CUR );
    fseek( input, 4*( inemdata->dims[0]*inemdata->dims[1]*( nr[2] - 1 ) ), SEEK_CUR );
    ilaufx = 0;
    ilaufz = 0;
    count = 0;
    fseek_merker = 0;
    s1 = 0;s2 = 0;s3 = 0;
    xy_dims = area_d[0] * area_d[1];
    for ( lauf = nr[2];lauf <= nr[2] + area[2];lauf++ ) {
      for ( laufy = nr[1];laufy <= nr[1] + area[1];laufy++ ) {
        count = fread ( &inemdata->floatdata[( ilaufz*xy_dims )+ilaufx*area_d[0]], 4, area_d[0], input );
        fseek( input, ( long )( sizeof( float )*( inemdata->dims[0] - area[0] - 1 ) ), SEEK_CUR );
        fseek_merker = fseek_merker + inemdata->dims[0];
        ilaufx = ilaufx + 1;
      }
      ilaufz = ilaufz + 1;
      ilaufx = 0;
      size_area = 0;
      fseek( input, 4*( inemdata->dims[1]*inemdata->dims[0] - fseek_merker ), SEEK_CUR );
      fseek_merker = 0;
    }
    /*  fread (inemdata->floatdata, 4, size, input); */
  }
  fclose ( input );

}

void write_em ( char *outfile, struct em_file *outemdata ) {
//     int i;
  long size;
  FILE *output = 0;

  if ( ( output = fopen ( outfile, "w" ) ) == 0 ) {
    printf( "could not open %s\n", outfile );
    exit ( 23 );
  }


  fwrite ( outemdata->magic, 1, 1, output );
  fwrite ( outemdata->dummya, 1, 2, output );
  fwrite ( outemdata->type, 1, 1, output );
  fwrite ( outemdata->dims, 4, 3, output );
  fwrite ( outemdata->comment, 1, 80, output );
  fwrite ( outemdata->emdata, 4, 40, output );
  fwrite ( outemdata->dummyb, 1, 256, output );


  size = outemdata->dims[0] * outemdata->dims[1] * outemdata->dims[2];

  switch ( outemdata->type[0] ) {
  case 1: fwrite ( outemdata->bytedata, 1, size, output );
  case 2: fwrite ( outemdata->intdata, 2, size, output );
  case 5: fwrite ( outemdata->floatdata, 4, size, output );
  }

  fclose ( output );
}

void write_em_subregion ( char *outfile, struct em_file *outemdata, int *nr, int *area ) {
  FILE *output = 0;
//   long size;
  int lauf, ilaufx, ilaufz, laufy ;
  long int s1, s2, s3;
  int area_d[3];
  long int xy_dims, fseek_merker;
//   unsigned int laufin;
  int size_area;
  int count;

  if ( ( output = fopen ( outfile, "r+b" ) ) == 0 ) {
    printf( "Could not open file for reading in write_em_subregion.\n" ); exit( 24 );
  }
  fread ( outemdata->magic, 1, 1, output );
  fread ( outemdata->dummya, 1, 2, output );
  fread ( outemdata->type, 1, 1, output );
  fread ( outemdata->dims, 4, 3, output );
  fread ( outemdata->comment, 1, 80, output );
  fread ( outemdata->emdata, 4, 40, output );
  fread ( outemdata->dummyb, 1, 256, output );
  if ( ( nr[0] + area[0] ) > outemdata->dims[0] + 1 || nr[1] + area[1] > outemdata->dims[1] + 1 || nr[2] + area[2] > outemdata->dims[2] + 1 )
    {printf( "Subregion dimensions plus offset larger than volume dimensions." ); return;}

  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;

  fseek( output, 4*( nr[0] - 1 ), SEEK_CUR );
  fseek( output, 4*( outemdata->dims[0]*( nr[1] - 1 ) ), SEEK_CUR );
  fseek( output, 4*( outemdata->dims[0]*outemdata->dims[1]*( nr[2] - 1 ) ), SEEK_CUR );
  ilaufx = 0;
  ilaufz = 0;
  count = 0;
  fseek_merker = 0;
  s1 = 0;s2 = 0;s3 = 0;
  xy_dims = area_d[0] * area_d[1];
  for ( lauf = nr[2];lauf <= nr[2] + area[2];lauf++ ) {
    for ( laufy = nr[1];laufy <= nr[1] + area[1];laufy++ ) {
      count = fwrite ( &outemdata->floatdata[( ilaufz*xy_dims )+ilaufx*area_d[0]], 4, area_d[0], output );
      fseek( output, ( long )( sizeof( float )*( outemdata->dims[0] - area[0] - 1 ) ), SEEK_CUR );
      fseek_merker = fseek_merker + outemdata->dims[0];
      ilaufx = ilaufx + 1;
    }
    ilaufz = ilaufz + 1;
    ilaufx = 0;
    size_area = 0;
    fseek( output, 4*( outemdata->dims[1]*outemdata->dims[0] - fseek_merker ), SEEK_CUR );
    fseek_merker = 0;

  }
  fclose ( output );
}

void write_em_subsubregion ( char *outfile, struct em_file *outemdata, int *nr, int *area, int *offset, int *dimarray ) {
  FILE *output = 0;
//   long size;
  int lauf, ilaufx, ilaufz, laufy ;
  long int s1, s2, s3;
  int area_d[3], ioffset[3];
  long int xy_dims, fseek_merker,/*xy_offset,*/xy_dimarray, initial_offset;
//   unsigned int laufin;
  int size_area;
  int count;
  /* outfile      filename
     outemdata    array to write out
     nr           start position
     area
     offset       offset in array - attention: offset 0 means NO offset, i.e. start at coord (1,1,1) (in Matlab convention!)
     dimarray     dimension of array
  */

  if ( ( output = fopen ( outfile, "r+b" ) ) == 0 ) {
    printf( "Could not open file for reading in write_em_subsubregion.\n" );  exit( 25 );
  }
  fread ( outemdata->magic, 1, 1, output );
  fread ( outemdata->dummya, 1, 2, output );
  fread ( outemdata->type, 1, 1, output );
  fread ( outemdata->dims, 4, 3, output );
  fread ( outemdata->comment, 1, 80, output );
  fread ( outemdata->emdata, 4, 40, output );
  fread ( outemdata->dummyb, 1, 256, output );
  if ( ( nr[0] + area[0] ) > outemdata->dims[0] + 1 || nr[1] + area[1] > outemdata->dims[1] + 1 || nr[2] + area[2] > outemdata->dims[2] + 1 )
    {printf( "Subregion dimensions plus offset larger than volume dimensions." ); return;}

  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;
  fseek( output, 4*( nr[0] - 1 ), SEEK_CUR );
  fseek( output, 4*( outemdata->dims[0]*( nr[1] - 1 ) ), SEEK_CUR );
  fseek( output, 4*( outemdata->dims[0]*outemdata->dims[1]*( nr[2] - 1 ) ), SEEK_CUR );
  ilaufx = 0;
  ilaufz = 0;
  count = 0;
  fseek_merker = 0;
  s1 = 0;s2 = 0;s3 = 0;
  ioffset[0] = offset[0];
  ioffset[1] = offset[1];
  ioffset[2] = offset[2];
  xy_dims = area_d[0] * area_d[1];
  xy_dimarray = dimarray[0] * dimarray[1];
  initial_offset = xy_dimarray * ioffset[2] + ioffset[1] * dimarray[1]; /* dimxy*z_start + dimx*y_start ; bug? dimarray[1] -> dimarray[2] FF*/
  /*printf("initial offset: %i , ioffset[0]: %i , ioffset[1]: %i , ioffset[2]: %i , xy_dimarray : %i \n",
    initial_offset,ioffset[0],ioffset[1],ioffset[2],xy_dimarray);fflush(stdout);*/
  for ( lauf = nr[2];lauf <= nr[2] + area[2];lauf++ ) {
    for ( laufy = nr[1];laufy <= nr[1] + area[1];laufy++ ) {
      count = fwrite ( &outemdata->floatdata[( ilaufz*xy_dimarray )+ilaufx*dimarray[1] + initial_offset + ioffset[0]], 4, area_d[0], output ); /*bug: ilauf feangt bei 1 an! FF */
      fseek( output, ( long )( sizeof( float )*( outemdata->dims[0] - area[0] - 1 ) ), SEEK_CUR );
      fseek_merker = fseek_merker + outemdata->dims[0];
      ilaufx = ilaufx + 1;
    }
    ilaufz = ilaufz + 1;
    ilaufx = 0;
    size_area = 0;
    fseek( output, 4*( outemdata->dims[1]*outemdata->dims[0] - fseek_merker ), SEEK_CUR );
    fseek_merker = 0;
  }
  fclose ( output );
}



void create_em ( char *outfile, int *nr ) {
  unsigned char magic[1];
//     int f_magic;
  char dummya[2];
  unsigned char type[1];
  int dims[3];
  char comment[80];
  int emdata[40];
  char dummyb[256];
  float *floatdata;
  FILE *output = 0;
//     long size;
  int lauf;

  if ( ( output = fopen ( outfile, "wb" ) ) == 0 ) {
    printf( "Could not create file in create_em.\n" ); exit( 26 );
  }
  magic[0] = 6; /* for PC */
  fwrite ( magic, 1, 1, output );
  dummya[0] = 0;
  dummya[1] = 0;
  fwrite ( dummya, 1, 2, output );
  type[0] = 5; /* for float */
  fwrite ( type, 1, 1, output );
  dims[0] = nr[0];
  dims[1] = nr[1];
  dims[2] = nr[2];
  fwrite ( dims, 4, 3, output );
  for ( lauf = 0;lauf < 80;lauf++ ){comment[lauf] = 0.0;}
  fwrite ( comment, 1, 80, output );
  for ( lauf = 0;lauf < 40;lauf++ ){emdata[lauf] = 0.0;}
  fwrite ( emdata, 4, 40, output );
  for ( lauf = 0;lauf < 256;lauf++ ){dummyb[lauf] = 0.0;}
  fwrite ( dummyb, 1, 256, output );
  if ( ( floatdata = ( float * )( malloc( dims[1] * dims[0] * sizeof( float ) ) ) ) == 0 )
    {printf( "Memory allocation problem in tom_emwritec.\n" );exit( 27 );}
  for ( lauf = 0;lauf < dims[1]*dims[0];lauf++ ){floatdata[lauf] = 0.0;};
  for ( lauf = 0;lauf < dims[2];lauf++ ){fwrite( &floatdata[0], sizeof( float ), dims[0]*dims[1], output );};
  fflush( output );
  fclose( output );
  free( floatdata );
  return;
}
