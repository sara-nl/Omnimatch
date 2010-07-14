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

#include "em.h"

#define EM_HEADER_LEN 512

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

enum MachineCode {OS9=0, VAX=1, CONVEX=2, SGI=3, SUN=4, MAC=5, PC=6};
enum TypeCode {BYTE=1, SHORT=2, LONG=4, FLOAT=5, COMPLEX=8, DOUBLE=9};


FILE *em_header_read ( const char *filename, struct gen_data *emdata, const char * mode) {
  FILE *file = 0;
  if ( ( file = fopen ( filename, mode ) ) == 0 ) {
    printf( "could not open %s in em_header_read\n", filename );
    exit ( 1 );
  }
  fread ( emdata->header.emHeader.magic, 1, 1, file );
  fread ( emdata->header.emHeader.dummya, 1, 2, file );
  fread ( emdata->header.emHeader.type, 1, 1, file );
  fread ( emdata->header.emHeader.dims, 4, 3, file );
  fread ( emdata->header.emHeader.comment, 1, 80, file );
  fread ( emdata->header.emHeader.emdata, 4, 40, file );
  fread ( emdata->header.emHeader.dummyb, 1, 256, file );
  return file;
};

int isEM( const char *filename) {
  FILE *input = 0;
  struct gen_data data;
  long position = 0;
  input = em_header_read(filename, &data, "r");
  fseek(input, 0L, SEEK_END);
  position = ftell( input );
  fclose ( input );
  switch (data.header.emHeader.type[0]) {
  case BYTE: return position == (data.header.emHeader.dims[0] * data.header.emHeader.dims[1] * data.header.emHeader.dims[2] + EM_HEADER_LEN);
  case SHORT: return position == ((int)sizeof(int) * data.header.emHeader.dims[0] * data.header.emHeader.dims[1] * data.header.emHeader.dims[2] + EM_HEADER_LEN);
  case FLOAT: return position == ((int)sizeof(float) * data.header.emHeader.dims[0] * data.header.emHeader.dims[1] * data.header.emHeader.dims[2] + EM_HEADER_LEN);
  default: return 0;
  }
}

void read_em_header ( const char *infile, struct gen_data *inemdata) {
  FILE *input = 0;
  input = em_header_read(infile, inemdata, "r");
  fclose ( input );
}

void read_em ( const char *infile, struct gen_data *inemdata ) {
  FILE *input = 0;
  long size;
  input = em_header_read(infile, inemdata, "r");
  size = inemdata->header.emHeader.dims[0] * inemdata->header.emHeader.dims[1] * inemdata->header.emHeader.dims[2];
  switch ( inemdata->header.emHeader.type[0] ) {
  case BYTE: if ( ( inemdata->bdata = ( unsigned char * )( malloc( size ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 1 ); }
    break;
  case SHORT: if ( ( inemdata->idata = ( int * )( malloc( size * 2 ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 1 ); }
    break;
  case FLOAT: if ( ( inemdata->fdata = ( float * )( malloc( size * 4 ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 1 ); }
    break;
  }
  switch ( inemdata->header.emHeader.type[0] ) {
  case BYTE: fread ( inemdata->bdata, 1, size, input );
  case SHORT: fread ( inemdata->idata, 2, size, input );
  case FLOAT: fread ( inemdata->fdata, 4, size, input );
  }
  fclose ( input );
};

void read_em_subregion ( const char *infile, struct gen_data *inemdata, int *nr, int  *area ) {
  FILE *input = 0;
  long size;
  int lauf, ilaufx, ilaufz, laufy ;
  int area_d[3];
  long int xy_dims, fseek_merker;
  int count;
  input = em_header_read(infile, inemdata, "r");
  size = ( nr[0] + area[0] ) * ( nr[1] + area[1] ) * ( nr[2] + area[2] );
  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;
  if (inemdata->header.emHeader.type[0] == FLOAT) {
    fseek( input, 4*(inemdata->header.emHeader.dims[0]*inemdata->header.emHeader.dims[1]*(nr[2] - 1) + inemdata->header.emHeader.dims[0]*(nr[1] - 1) + nr[0] - 1), SEEK_CUR );
    ilaufx = 0;
    ilaufz = 0;
    count = 0;
    fseek_merker = 0;
    xy_dims = area_d[0] * area_d[1];
    for ( lauf = nr[2]; lauf <= nr[2] + area[2]; lauf++ ) {
      for ( laufy = nr[1]; laufy <= nr[1] + area[1]; laufy++ ) {
        count = fread ( &inemdata->fdata[( ilaufz*xy_dims )+ilaufx*area_d[0]], 4, area_d[0], input );
        fseek( input, ( long )( sizeof( float )*( inemdata->header.emHeader.dims[0] - area[0] - 1 ) ), SEEK_CUR );
        fseek_merker = fseek_merker + inemdata->header.emHeader.dims[0];
        ilaufx++;
      }
      ilaufz++;
      ilaufx = 0;
      fseek( input, 4*( inemdata->header.emHeader.dims[1]*inemdata->header.emHeader.dims[0] - fseek_merker ), SEEK_CUR );
      fseek_merker = 0;
    }
  }
  fclose ( input );
}

void write_em ( const char *outfile, struct gen_data *outemdata ) {
  long size;
  FILE *output = 0;
  output = em_header_read(outfile, outemdata, "w");
  size = outemdata->header.emHeader.dims[0] * outemdata->header.emHeader.dims[1] * outemdata->header.emHeader.dims[2];
  switch ( outemdata->header.emHeader.type[0] ) {
  case BYTE: fwrite ( outemdata->bdata, 1, size, output );
  case SHORT: fwrite ( outemdata->idata, 2, size, output );
  case FLOAT: fwrite ( outemdata->fdata, 4, size, output );
  }
  fclose ( output );
}

void write_em_subregion ( const char *outfile, struct gen_data *outemdata, int *nr, int *area ) {
  FILE *output = 0;
  int lauf, ilaufx, ilaufz, laufy ;
  int area_d[3];
  long int xy_dims, fseek_merker;
  int count;
  output = em_header_read(outfile, outemdata, "r+b");
  if ( ( nr[0] + area[0] ) > outemdata->header.emHeader.dims[0] + 1 || nr[1] + area[1] > outemdata->header.emHeader.dims[1] + 1 || nr[2] + area[2] > outemdata->header.emHeader.dims[2] + 1 ) 
    {printf( "Subregion dimensions plus offset larger than volume dimensions." ); return;}
  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;
  fseek( output, 4*(outemdata->header.emHeader.dims[0]*outemdata->header.emHeader.dims[1]*(nr[2] - 1) + outemdata->header.emHeader.dims[0]*(nr[1] - 1) + nr[0] - 1), SEEK_CUR );
  ilaufx = 0;
  ilaufz = 0;
  count = 0;
  fseek_merker = 0;
  xy_dims = area_d[0] * area_d[1];
  for ( lauf = nr[2]; lauf <= nr[2] + area[2]; lauf++ ) {
    for ( laufy = nr[1]; laufy <= nr[1] + area[1]; laufy++ ) {
      count = fwrite ( &outemdata->fdata[( ilaufz*xy_dims )+ilaufx*area_d[0]], 4, area_d[0], output );
      fseek( output, ( long )( sizeof( float )*( outemdata->header.emHeader.dims[0] - area[0] - 1 ) ), SEEK_CUR );
      fseek_merker = fseek_merker + outemdata->header.emHeader.dims[0];
      ilaufx++;
    }
    ilaufz++;
    ilaufx = 0;
     fseek( output, 4*( outemdata->header.emHeader.dims[1]*outemdata->header.emHeader.dims[0] - fseek_merker ), SEEK_CUR );
    fseek_merker = 0;
  }
  fclose ( output );
}

/* outfile      filename
   outemdata    array to write out
   nr           start position
   area
   offset       offset in array - attention: offset 0 means NO offset, i.e. start at coord (1,1,1) (in Matlab convention!)
   dimarray     dimension of array
*/
void write_em_subsubregion ( const char *outfile, struct gen_data *outemdata, int *nr, int *area, int *offset, int *dimarray ) {
  FILE *output = 0;
  int lauf, ilaufx, ilaufz, laufy ;
  int area_d[3];
  long int fseek_merker, xy_dimarray, initial_offset;
  int count;
  output = em_header_read(outfile, outemdata, "r+b");
  if ( ( nr[0] + area[0] ) > outemdata->header.emHeader.dims[0] + 1 || nr[1] + area[1] > outemdata->header.emHeader.dims[1] + 1 || nr[2] + area[2] > outemdata->header.emHeader.dims[2] + 1 )
    {printf( "Subregion dimensions plus offset larger than volume dimensions." ); return;}
  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;
  fseek( output, 4*(outemdata->header.emHeader.dims[0]*outemdata->header.emHeader.dims[1]*(nr[2] - 1) + outemdata->header.emHeader.dims[0]*(nr[1] - 1) + nr[0] - 1), SEEK_CUR );
  ilaufx = 0;
  ilaufz = 0;
  count = 0;
  fseek_merker = 0;
  xy_dimarray = dimarray[0] * dimarray[1];
  initial_offset = xy_dimarray * offset[2] + offset[1] * dimarray[1]; /* dimxy*z_start + dimx*y_start ; bug? dimarray[1] -> dimarray[2] FF*/
  /*printf("initial offset: %i , ioffset[0]: %i , ioffset[1]: %i , ioffset[2]: %i , xy_dimarray : %i \n",
    initial_offset,ioffset[0],ioffset[1],ioffset[2],xy_dimarray);fflush(stdout);*/
  for ( lauf = nr[2]; lauf <= nr[2] + area[2]; lauf++ ) {
    for ( laufy = nr[1]; laufy <= nr[1] + area[1]; laufy++ ) {
      count = fwrite ( &outemdata->fdata[( ilaufz*xy_dimarray )+ilaufx*dimarray[1] + initial_offset + offset[0]], 4, area_d[0], output ); /*bug: ilauf feangt bei 1 an! FF */
      fseek( output, ( long )( sizeof( float )*( outemdata->header.emHeader.dims[0] - area[0] - 1 ) ), SEEK_CUR );
      fseek_merker = fseek_merker + outemdata->header.emHeader.dims[0];
      ilaufx++;
    }
    ilaufz++;
    ilaufx = 0;
    fseek( output, 4*( outemdata->header.emHeader.dims[1]*outemdata->header.emHeader.dims[0] - fseek_merker ), SEEK_CUR );
    fseek_merker = 0;
  }
  fclose ( output );
}



void create_em ( const char *outfile, int *nr ) {
  unsigned char magic[1];
  char dummya[2];
  unsigned char type[1];
  int dims[3];
  char comment[80];
  int emdata[40];
  char dummyb[256];
  float *floatdata;
  FILE *output = 0;
  int lauf;
  if ( ( output = fopen ( outfile, "wb" ) ) == 0 ) {
    printf( "Could not create file in create_em.\n" ); exit( 1 );
  }
  magic[0] = PC; /* for PC */
  fwrite ( magic, 1, 1, output );
  dummya[0] = 0;
  dummya[1] = 0;
  fwrite ( dummya, 1, 2, output );
  type[0] = FLOAT; /* for float */
  fwrite ( type, 1, 1, output );
  dims[0] = nr[0];
  dims[1] = nr[1];
  dims[2] = nr[2];
  fwrite ( dims, 4, 3, output );
  for ( lauf = 0; lauf < 80; lauf++ ) 
    comment[lauf] = 0.0;
  fwrite ( comment, 1, 80, output );
  for ( lauf = 0; lauf < 40; lauf++ ) 
    emdata[lauf] = 0.0;
  fwrite ( emdata, 4, 40, output );
  for ( lauf = 0; lauf < 256; lauf++ )
    dummyb[lauf] = 0.0;
  fwrite ( dummyb, 1, 256, output );
  if ( ( floatdata = ( float * )( malloc( dims[1] * dims[0] * sizeof( float ) ) ) ) == 0 )
    {printf( "Memory allocation problem in tom_emwritec.\n" );exit( 1 );}
  for ( lauf = 0; lauf < dims[1]*dims[0]; lauf++ )
    floatdata[lauf] = 0.0;
  for ( lauf = 0; lauf < dims[2]; lauf++ ) 
    fwrite( &floatdata[0], sizeof( float ), dims[0]*dims[1], output );
  fflush( output );
  fclose( output );
  free( floatdata );
}
