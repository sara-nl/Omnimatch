#include <rfftw.h>
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

void read_em (char *infile, struct em_file *inemdata);
void read_em_header (char *infile, struct em_file *inemdata);
void read_em_subregion (char *infile, struct em_file *inemdata, int *nr, int *area);
void write_em (char *outfile, struct em_file *outemdata);
void write_em_subregion (char *infile, struct em_file *inemdata, int *nr, int *area);
void write_em_subsubregion (char *outfile, struct em_file *outemdata, int *nr, int *area,int *offset,int *dimarray);
// pastes.c        by S.N. modified by F.F.
void pastes(float *I, float *O,int Ox_min, int Oy_min, int Oz_min, int Ox_max, int Oy_max, int Oz_max,int Vx_max);
// energizer.c     by F.F.
float energizer (int, int, int, float *, float *, float *, rfftwnd_plan, rfftwnd_plan);
int count_voxel (int, float *, float);
void cross(float *volinout, int Vx_max);
void shift(fftw_complex *vol_four, float dx, float dy, float dz, int Nx, int Ny, int Nz);
float sumvoxel (int, int, int, float *);/* sums voxels > eps of inputarray  */
float prepref(int Rx_min, int Rx_max, float eps, float *Rot_tmpl, float *inputdata, float *mask, rfftwnd_plan, rfftwnd_plan);
/* prepare reference:  */
/* sort4fftw  */
void sort4fftw(fftw_real *ext_array, float *floatdata, int Nx, int Ny, int Nz);
void sortback4fftw(fftw_real *ext_array, float *floatdata, int Nx, int Ny, int Nz);
/* four_filter - filters in Fourier space */
void lowpass(fftw_complex *vol_four, float R, float smooth, int Nx, int Ny, int Nz, fftw_real scale);
void bandpass(fftw_complex *vol_four, float R_down, float R_up, float smooth, int Nx, int Ny, int Nz, fftw_real scale);
void correl(fftw_complex *vol_four, fftw_complex *ref_four, int Nx, int Ny, int Nz, fftw_real scale);
void convolve(fftw_complex *vol_four, fftw_complex *ref_four, int Nx, int Ny, int Nz, fftw_real scale);
/* correlation is written to ref_four */
/* real_utils - real space utilities */
void symref(float *volume, int nfold, int Nx, int Ny, int Nz);
float variance(float *volume, int Nx, int Ny, int Nz);
void norm(float *volume, int Nx, int Ny, int Nz);
void limit(float *volume, int Nx, int Ny, int Nz, float low, float hi);
void limitz(float *volume, int Nx, int Ny, int Nz, float low, float hi);
void pointrotate(float x, float y, float z, float phi, float psi, float theta);

