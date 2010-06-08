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

//#define USE_GPUS
#include <sarafft.h>

// tom_rotate3d.c  by S.N.
void tom_rotate3d( float *O, float *I, float Phi, float Psi, float The, int Ox_max, int Oy_max, int Oz_max);
void rotate3d( float *O, float *I, float Phi, float Psi, float The, int Ox_max, int Oy_max, int Oz_max);
// emfile.c        by S.N.
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
void create_em (char *outfile, int *nr);
void read_em (char *infile, struct em_file *inemdata);
void read_em_header (char *infile, struct em_file *inemdata);
void read_em_subregion (char *infile, struct em_file *inemdata, int *nr, int *area);
void write_em (char *outfile, struct em_file *outemdata);
void write_em_subregion (char *infile, struct em_file *inemdata, int *nr, int *area);
void write_em_subsubregion (char *outfile, struct em_file *outemdata, int *nr, int *area,int *offset,int *dimarray);
// pastes.c        by S.N. modified by F.F.
void pastes(float *I, float *O,int Ox_min, int Oy_min, int Oz_min, int Ox_max, int Oy_max, int Oz_max,int Vx_max);
// energizer.c     by F.F.
// PieterB: The original declaration had an int as the third parameter, but this was not consistent with the implementation!
//float energizer (int, int, int, float *, float *, float *, sararfftnd_plan, sararfftnd_plan);
float energizer (int, int, float, float *, float *, float *, sararfftnd_plan, sararfftnd_plan);
int count_voxel (int, float *, float);
void cross(float *volinout, int Vx_max);
void shift(sarafft_complex *vol_four, float dx, float dy, float dz, int Nx, int Ny, int Nz);
float sumvoxel (int, int, int, float *);/* sums voxels > eps of inputarray  */
float prepref(int Rx_min, int Rx_max, float eps, float *Rot_tmpl, float *inputdata, float *mask, sararfftnd_plan, sararfftnd_plan);
/* prepare reference:  */
/* sort4fftw  */
void sort4fftw(sarafft_real *ext_array, float *floatdata, int Nx, int Ny, int Nz);
void sortback4fftw(sarafft_real *ext_array, float *floatdata, int Nx, int Ny, int Nz);
/* four_filter - filters in Fourier space */
void lowpass(sarafft_complex *vol_four, float R, float smooth, int Nx, int Ny, int Nz, sarafft_real scale);
void bandpass(sarafft_complex *vol_four, float R_down, float R_up, float smooth, int Nx, int Ny, int Nz, sarafft_real scale);
void correl(sarafft_complex *vol_four, sarafft_complex *ref_four, int Nx, int Ny, int Nz, sarafft_real scale);
void convolve(sarafft_complex *vol_four, sarafft_complex *ref_four, int Nx, int Ny, int Nz, sarafft_real scale);
/* correlation is written to ref_four */
/* real_utils - real space utilities */
void symref(float *volume, int nfold, int Nx, int Ny, int Nz);
float variance(float *volume, int Nx, int Ny, int Nz);
void norm(float *volume, int Nx, int Ny, int Nz);
void limit(float *volume, int Nx, int Ny, int Nz, float low, float hi);
void limitz(float *volume, int Nx, int Ny, int Nz, float low, float hi);
/*void pointrotate(float x, float y, float z, float phi, float psi, float theta);*/

