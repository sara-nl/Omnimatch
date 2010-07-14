#ifndef DATATYPES_H
#define DATATYPES_H

/* 
Structure of EM-data files: (from Frangakis group website)
This module handles EM files (Hegerl, R. (1996). The EM Program Package: A Platform for Image Processing in Biological Electron Microscopy, J. Struct. Biol.  116: 30-34), 
with the default extension .em. The file can be in little endian or big endian byte ordering with a 512 byte header. The structure of an em file is as follows:
- Byte1: Machine coding 	Machine 	Value
  	OS-9 	0
  	VAX 	1
  	CONVEX 	2
  	SGI 	3
  	SUN 	4 (not supported)
  	MAC 	5
  	PC 	6
- Byte 2: General purpose 	On OS-9 	old version 0, new version 1
- Byte 3: Not used in EM-format 	if 1, the header is abandonned
- Byte 4: Data type coding 	Image type 	No. of byte 	Value 	 
  	byte 	1 	1 	 
  	short 	2 	2 	 
  	long int 	4 	4 	 
  	float 	4 	5 	 
  	complex 	8 	8 	 
  	double 	8 	9 	 
- Three long integers (3x4 bytes) are image size in x, y, z dimension
- 80 characters as comment
- 40 long integers (4x40 bytes) are user defined parameters
- 256 bytes with userdata, first 20 chars username, 8 chars date (i.e. 03/02/03)
- raw data following with the x variable as the fastest dimensions, then y and z
*/

struct em_header {
  unsigned char magic[1];
  char dummya[2];
  unsigned char type[1];
  int dims[3];
  char comment[80];
  int emdata[40];
  char dummyb[256];
};

/*
Structure of MRC-data files:
MRC Header has a length of 1024 bytes
 SIZE  DATA    NAME    DESCRIPTION
   4   int     NX      number of Columns    (fastest changing in map)
   4   int     NY      number of Rows
   4   int     NZ      number of Sections   (slowest changing in map)
   4   int     MODE    Types of pixel in image
                       0 = Image     unsigned bytes
                       1 = Image     signed short integer (16 bits)
                       2 = Image     float
                       3 = Complex   short*2
                       4 = Complex   float*2
   4   int     NXSTART Number of first COLUMN  in map (Default = 0)
   4   int     NYSTART Number of first ROW     in map      "
   4   int     NZSTART Number of first SECTION in map      "
   4   int     MX      Number of intervals along X
   4   int     MY      Number of intervals along Y
   4   int     MZ      Number of intervals along Z
   4   float   Xlen    Cell Dimensions (Angstroms)
   4   float   Ylen                 "
   4   float   Zlen                 "
   4   float   ALPHA   Cell Angles (Degrees)
   4   float   BETA                 "
   4   float   GAMMA                "
   4   int     MAPC    Which axis corresponds to Columns  (1,2,3 for X,Y,Z)
   4   int     MAPR    Which axis corresponds to Rows     (1,2,3 for X,Y,Z)
   4   int     MAPS    Which axis corresponds to Sections (1,2,3 for X,Y,Z)
   4   float   AMIN    Minimum density value
   4   float   AMAX    Maximum density value
   4   float   AMEAN   Mean    density value    (Average)
   2   short   ISPG    Space group number       (0 for images)
   2   short   NSYMBT  Number of bytes used for storing symmetry operators
   4   int     NEXT    Number of bytes in extended header
   2   short   CREATID Creator ID
   30    -     EXTRA   Not used. All set to zero by default
   2   short   NINT    Number of integer per section
   2   short   NREAL   Number of reals per section
   28    -     EXTRA2  Not used. All set to zero by default
   2   short   IDTYPE  0=mono, 1=tilt, 2=tilts, 3=lina, 4=lins
   2   short   LENS
   2   short   ND1
   2   short   ND2
   2   short   VD1
   2   short   VD2
   24  float   TILTANGLES
   4   float   XORIGIN X origin
   4   float   YORIGIN Y origin
   4   float   ZORIGIN Z origin
   4   char    CMAP    Contains "MAP "
   4   char    STAMP
   4   float   RMS
   4   int     NLABL   Number of labels being used
   800 char    10 labels of 80 character

HEADER EXTENSION [OPTIONAL]
 SIZE  DATA    NAME    DESCRIPTION
   4   float   a_tilt  Tilt angle
   4   float   b_tilt  Tilt angle
   4   float   x_stage ?
   4   float   y_stage ?
   4   float   z_stage ?
   4   float   x_shift ?
   4   float   y_shift ?
   4   float   defocus ?
   4   float   exp_time Exposure time
   4   float   mean_int ?
   4   float   tiltaxis Tilt axis
   4   float   pixelsize Object pixel size
   4   float   magnification Magnification

*/

/* enum {MRC_BYTE=0, MRC_SHORT, MRC_FLOAT, MRC_COMPLEX1, MRC_COMPLEX2}; */
/* enum TypeCode {BYTE=1, SHORT=2, LONG=4, FLOAT=5, COMPLEX=8, DOUBLE=9}; */

struct mrc_header {
  int   nx[1];         
  int   ny[1];        
  int   nz[1];         
  int   mode[1];       
  int   nxstart[1];    
  int   nystart[1];
  int   nzstart[1];
  int   mx[1];        
  int   my[1];
  int   mz[1];
  float   xlen[1];
  float   ylen[1];       
  float   zlen[1];
  float   alpha[1];      
  float   beta[1];
  float   gamma[1];
  int   mapc[1];       
  int   mapr[1];       
  int   maps[1];       
  float   amin[1];
  float   amax[1];
  float   amean[1];
  short int   ispg[1];      
  short int   nsymbt[1];     
  int   next[1];
  short int   creatid[1];  
  char    extra[30];
  short int   nint[1];
  short int   nreal[1];
  char    extra2[28];
  short int   idtype[1];
  short int   lens[1];
  short int   nd1[1];     
  short int   nd2[1];
  short int   vd1[1];
  short int   vd2[1];
  float   tiltangles[6];  
  float   xorg[1];
  float   yorg[1];
  float   zorg[1];
  char    cmap[4];
  char    stamp[4];
  float   rms[1];
  int     nlabl[1];
  char    labels[800];
};

union gen_header {
  struct em_header emHeader;
  struct mrc_header mrcHeader;
};

struct gen_data {
  union gen_header header;
  unsigned char *bdata;
  int *idata;
  float *fdata;  
};


#endif
