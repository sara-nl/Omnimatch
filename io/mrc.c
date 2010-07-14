
#include <stdio.h>
#include <stdlib.h>

#include "mrc.h"

#define MRC_HEADER_LEN 1024

enum {MRC_BYTE=0, MRC_SHORT, MRC_FLOAT, MRC_COMPLEX1, MRC_COMPLEX2};


FILE *mrc_header_read ( const char *filename, struct gen_data *mrcdata, const char *mode) {
  FILE *file = 0;
  if ( ( file = fopen ( filename, mode ) ) == 0 ) {
    printf( "could not open %s in mrc_header_read\n", filename );
    exit ( 1 );
  }

  /* Read dimensions */
  fread ( mrcdata->header.mrcHeader.nx, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.ny, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.nz, 4, 1, file );
  /* Read data type */
  fread ( mrcdata->header.mrcHeader.mode, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.nxstart, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.nystart, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.nzstart, 4, 1, file );
  
  fread ( mrcdata->header.mrcHeader.mx, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.my, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.mz, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.xlen, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.ylen, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.zlen, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.alpha, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.beta, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.gamma, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.mapc, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.mapr, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.maps, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.amin, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.amax, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.amean, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.ispg, 2, 1, file );
  fread ( mrcdata->header.mrcHeader.nsymbt, 2, 1, file );

  fread ( mrcdata->header.mrcHeader.next, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.creatid, 2, 1, file );

  fread ( mrcdata->header.mrcHeader.extra, 1, 30, file );

  fread ( mrcdata->header.mrcHeader.nint, 2, 1, file );
  fread ( mrcdata->header.mrcHeader.nreal, 2, 1, file );

  fread ( mrcdata->header.mrcHeader.extra2, 1, 28, file );

  fread ( mrcdata->header.mrcHeader.idtype, 2, 1, file );
  fread ( mrcdata->header.mrcHeader.lens, 2, 1, file );
  fread ( mrcdata->header.mrcHeader.nd1, 2, 1, file );
  fread ( mrcdata->header.mrcHeader.nd2, 2, 1, file );
  fread ( mrcdata->header.mrcHeader.vd1, 2, 1, file );
  fread ( mrcdata->header.mrcHeader.vd2, 2, 1, file );

  fread ( mrcdata->header.mrcHeader.tiltangles, 4, 6, file );

  fread ( mrcdata->header.mrcHeader.xorg, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.yorg, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.zorg, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.cmap, 1, 4, file );
  fread ( mrcdata->header.mrcHeader.stamp, 1, 4, file );
  fread ( mrcdata->header.mrcHeader.rms, 4, 1, file );

  fread ( mrcdata->header.mrcHeader.nlabl, 4, 1, file );
  fread ( mrcdata->header.mrcHeader.labels, 1, 800, file );

  return file;
};

int isMRC( const char *filename) {
  FILE *input = 0;
  struct gen_data data;
  long position = 0;
  input = mrc_header_read(filename, &data, "r");
  fseek(input, 0L, SEEK_END);
  position = ftell( input );
  fclose ( input );
  switch (data.header.mrcHeader.mode[0]) {
  case MRC_BYTE: return position == (data.header.mrcHeader.nx[0] * data.header.mrcHeader.ny[0] * data.header.mrcHeader.nz[0] + MRC_HEADER_LEN);
  case MRC_SHORT: return position == ((int)sizeof(int) * data.header.mrcHeader.nx[0] * data.header.mrcHeader.ny[0] * data.header.mrcHeader.nz[0] + MRC_HEADER_LEN);
  case MRC_FLOAT: return position == ((int)sizeof(float) * data.header.mrcHeader.nx[0] * data.header.mrcHeader.ny[0] * data.header.mrcHeader.nz[0] + MRC_HEADER_LEN);
  default: return 0;
  }
}

void read_mrc_header ( const char *infile, struct gen_data *mrcdata) {
  FILE *input = 0;
  input = mrc_header_read(infile, mrcdata, "r");
  fclose ( input );
}

void read_mrc ( const char *infile, struct gen_data *mrcdata ) {
  FILE *input = 0;
  long size;
  input = mrc_header_read(infile, mrcdata, "r");
  size = mrcdata->header.mrcHeader.nx[0] * mrcdata->header.mrcHeader.ny[0] * mrcdata->header.mrcHeader.nz[0];
  switch ( mrcdata->header.mrcHeader.mode[0] ) {
  case MRC_BYTE: if ( ( mrcdata->bdata = ( unsigned char * )( malloc( size ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 1 ); }
    break;
  case MRC_SHORT: if ( ( mrcdata->idata = ( int * )( malloc( size * 2 ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 1 ); }
    break;
  case MRC_FLOAT: if ( ( mrcdata->fdata = ( float * )( malloc( size * 4 ) ) ) == 0 )
      { printf ( "could not allocate memory" ); exit ( 1 ); }
    break;
  }
  switch ( mrcdata->header.mrcHeader.mode[0] ) {
  case MRC_BYTE: fread ( mrcdata->bdata, 1, size, input );
  case MRC_SHORT: fread ( mrcdata->idata, 2, size, input );
  case MRC_FLOAT: fread ( mrcdata->fdata, 4, size, input );
  }
  fclose ( input );
};

void read_mrc_subregion ( const char *infile, struct gen_data *mrcdata, int *nr, int  *area ) {
  FILE *input = 0;
  long size;
  int lauf, ilaufx, ilaufz, laufy ;
  int area_d[3];
  long int xy_dims, fseek_merker;
  int count;
  input = mrc_header_read(infile, mrcdata, "r");
  size = ( nr[0] + area[0] ) * ( nr[1] + area[1] ) * ( nr[2] + area[2] );
  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;
  if (mrcdata->header.mrcHeader.mode[0] == MRC_FLOAT) {
    fseek( input, 4*(mrcdata->header.mrcHeader.nx[0]*mrcdata->header.mrcHeader.ny[0]*(nr[2] - 1) + mrcdata->header.mrcHeader.nx[0]*(nr[1] - 1) + nr[0] - 1), SEEK_CUR );
    ilaufx = 0;
    ilaufz = 0;
    count = 0;
    fseek_merker = 0;
    xy_dims = area_d[0] * area_d[1];
    for ( lauf = nr[2]; lauf <= nr[2] + area[2]; lauf++ ) {
      for ( laufy = nr[1]; laufy <= nr[1] + area[1]; laufy++ ) {
        count = fread ( &mrcdata->fdata[( ilaufz*xy_dims )+ilaufx*area_d[0]], 4, area_d[0], input );
        fseek( input, ( long )( sizeof( float )*( mrcdata->header.mrcHeader.nx[0] - area[0] - 1 ) ), SEEK_CUR );
        fseek_merker = fseek_merker + mrcdata->header.mrcHeader.nx[0];
        ilaufx++;
      }
      ilaufz++;
      ilaufx = 0;
      fseek( input, 4*( mrcdata->header.mrcHeader.ny[0]*mrcdata->header.mrcHeader.nx[0] - fseek_merker ), SEEK_CUR );
      fseek_merker = 0;
    }
  }
  fclose ( input );
}

void write_mrc ( const char *outfile, struct gen_data *mrcdata ) {
  long size;
  FILE *output = 0;
  output = mrc_header_read(outfile, mrcdata, "w");
  size = mrcdata->header.mrcHeader.nx[0] * mrcdata->header.mrcHeader.ny[0] * mrcdata->header.mrcHeader.nz[0];
  switch ( mrcdata->header.mrcHeader.mode[0] ) {
  case MRC_BYTE: fwrite ( mrcdata->bdata, 1, size, output );
  case MRC_SHORT: fwrite ( mrcdata->idata, 2, size, output );
  case MRC_FLOAT: fwrite ( mrcdata->fdata, 4, size, output );
  }
  fclose ( output );
}

void write_mrc_subregion ( const char *outfile, struct gen_data *mrcdata, int *nr, int *area ) {
  FILE *output = 0;
  int lauf, ilaufx, ilaufz, laufy ;
  int area_d[3];
  long int xy_dims, fseek_merker;
  int count;
  output = mrc_header_read(outfile, mrcdata, "r+b");
  if ( ( nr[0] + area[0] ) > mrcdata->header.mrcHeader.nx[0] + 1 || nr[1] + area[1] > mrcdata->header.mrcHeader.ny[0] + 1 || nr[2] + area[2] > mrcdata->header.mrcHeader.nz[0] + 1 ) 
    {printf( "Subregion dimensions plus offset larger than volume dimensions." ); return;}
  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;
  fseek( output, 4*(mrcdata->header.mrcHeader.nx[0]*mrcdata->header.mrcHeader.ny[0]*(nr[2] - 1) + mrcdata->header.mrcHeader.nx[0]*(nr[1] - 1) + nr[0] - 1), SEEK_CUR );
  ilaufx = 0;
  ilaufz = 0;
  count = 0;
  fseek_merker = 0;
  xy_dims = area_d[0] * area_d[1];
  for ( lauf = nr[2]; lauf <= nr[2] + area[2]; lauf++ ) {
    for ( laufy = nr[1]; laufy <= nr[1] + area[1]; laufy++ ) {
      count = fwrite ( &mrcdata->fdata[( ilaufz*xy_dims )+ilaufx*area_d[0]], 4, area_d[0], output );
      fseek( output, ( long )( sizeof( float )*( mrcdata->header.mrcHeader.nx[0] - area[0] - 1 ) ), SEEK_CUR );
      fseek_merker = fseek_merker + mrcdata->header.mrcHeader.nx[0];
      ilaufx++;
    }
    ilaufz++;
    ilaufx = 0;
    fseek( output, 4*( mrcdata->header.mrcHeader.ny[0]*mrcdata->header.mrcHeader.nx[0] - fseek_merker ), SEEK_CUR );
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
void write_mrc_subsubregion ( const char *outfile, struct gen_data *mrcdata, int *nr, int *area, int *offset, int *dimarray ) {
  FILE *output = 0;
  int lauf, ilaufx, ilaufz, laufy ;
  int area_d[3];
  long int fseek_merker, xy_dimarray, initial_offset;
  int count;
  output = mrc_header_read(outfile, mrcdata, "r+b");
  if ( ( nr[0] + area[0] ) > mrcdata->header.mrcHeader.nx[0] + 1 || nr[1] + area[1] > mrcdata->header.mrcHeader.ny[0] + 1 || nr[2] + area[2] > mrcdata->header.mrcHeader.nz[0] + 1 )
    {printf( "Subregion dimensions plus offset larger than volume dimensions." ); return;}
  area_d[0] = area[0] + 1;
  area_d[1] = area[1] + 1;
  area_d[2] = area[2] + 1;
  fseek( output, 4*(mrcdata->header.mrcHeader.nx[0]*mrcdata->header.mrcHeader.ny[0]*(nr[2] - 1) + mrcdata->header.mrcHeader.nx[0]*(nr[1] - 1) + nr[0] - 1), SEEK_CUR );
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
      count = fwrite ( &mrcdata->fdata[( ilaufz*xy_dimarray )+ilaufx*dimarray[1] + initial_offset + offset[0]], 4, area_d[0], output ); /*bug: ilauf feangt bei 1 an! FF */
      fseek( output, ( long )( sizeof( float )*( mrcdata->header.mrcHeader.nx[0] - area[0] - 1 ) ), SEEK_CUR );
      fseek_merker = fseek_merker + mrcdata->header.mrcHeader.nx[0];
      ilaufx++;
    }
    ilaufz++;
    ilaufx = 0;
    fseek( output, 4*( mrcdata->header.mrcHeader.ny[0] * mrcdata->header.mrcHeader.nx[0] - fseek_merker ), SEEK_CUR );
    fseek_merker = 0;
  }
  fclose ( output );
}

void create_mrc ( const char *outfile, int *nr ) {
  float *floatdata;
  FILE *output = 0;
  struct mrc_header *header = 0;
  int counter;
  if ( ( output = fopen ( outfile, "wb" ) ) == 0 ) {
    printf( "Could not create file in create_mrc.\n" ); 
    exit( 1 );
  }

  header = (struct mrc_header*)calloc(1, sizeof(struct mrc_header));
  
  header->nx[0] = nr[0];
  header->ny[0] = nr[1];
  header->nz[0] = nr[2];

  fwrite ( header->nx, 4, 1, output );
  fwrite ( header->ny, 4, 1, output );
  fwrite ( header->nz, 4, 1, output );

  header->mode[0] = MRC_FLOAT;
  fwrite ( header->mode, 4, 1, output );

  fwrite ( header->nxstart, 4, 1, output );
  fwrite ( header->nystart, 4, 1, output );
  fwrite ( header->nzstart, 4, 1, output );

  fwrite ( header->mx, 4, 1, output );
  fwrite ( header->my, 4, 1, output );
  fwrite ( header->mz, 4, 1, output );

  fwrite ( header->xlen, 4, 1, output );
  fwrite ( header->ylen, 4, 1, output );
  fwrite ( header->zlen, 4, 1, output );

  fwrite ( header->alpha, 4, 1, output );
  fwrite ( header->beta, 4, 1, output );
  fwrite ( header->gamma, 4, 1, output );

  fwrite ( header->mapc, 4, 1, output );
  fwrite ( header->mapr, 4, 1, output );
  fwrite ( header->maps, 4, 1, output );

  fwrite ( header->amin, 4, 1, output );
  fwrite ( header->amax, 4, 1, output );
  fwrite ( header->amean, 4, 1, output );

  fwrite ( header->ispg, 2, 1, output );
  fwrite ( header->nsymbt, 2, 1, output );

  fwrite ( header->next, 4, 1, output );

  fwrite ( header->creatid, 2, 1, output );

  fwrite ( header->extra, 1, 30, output );

  fwrite ( header->nint, 2, 1, output );
  fwrite ( header->nreal, 2, 1, output );

  fwrite ( header->extra2, 1, 28, output );

  fwrite ( header->idtype, 2, 1, output );
  fwrite ( header->lens, 2, 1, output );

  fwrite ( header->nd1, 2, 1, output );
  fwrite ( header->nd2, 2, 1, output );

  fwrite ( header->vd1, 2, 1, output );
  fwrite ( header->vd2, 2, 1, output );

  fwrite ( header->tiltangles, 4, 6, output );

  fwrite ( header->xorg, 4, 1, output );
  fwrite ( header->yorg, 4, 1, output );
  fwrite ( header->zorg, 4, 1, output );

  fwrite ( header->cmap, 1, 4, output );
  fwrite ( header->stamp, 1, 4, output );

  fwrite ( header->rms, 4, 1, output );

  fwrite ( header->nlabl, 4, 1, output );

  fwrite ( header->labels, 1, 800, output );

  if ( ( floatdata = ( float * )( malloc( header->nx[0] * header->ny[0] * sizeof( float ) ) ) ) == 0 ) {
    printf( "Memory allocation problem in create_mrc.\n" );
    exit( 1 );
  }

  for ( counter = 0; counter < header->nx[0] * header->ny[0]; counter++ )
    floatdata[counter] = 0.0;
  for ( counter = 0; counter < header->nz[0]; counter++ ) 
    fwrite( &floatdata[0], sizeof( float ), header->nx[0] * header->ny[0], output );

  fflush( output );
  fclose( output );
  free( floatdata );
  free( header );
}

