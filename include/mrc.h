#ifndef MRC_H
#define MRC_H

#include "datatypes.h"

int  isMRC( const char *filename);
void create_mrc (const char *outfile, int *nr);
void read_mrc (const char *infile, struct gen_data *mrcdata);
void read_mrc_header (const char *infile, struct gen_data *mrcdata);
void read_mrc_subregion (const char *infile, struct gen_data *mrcdata, int *nr, int *area);
void write_mrc (const char *outfile, struct gen_data *mrcdata);
void write_mrc_subregion (const char *infile, struct gen_data *mrcdata, int *nr, int *area);
void write_mrc_subsubregion (const char *outfile, struct gen_data *mrcdata, int *nr, int *area,int *offset,int *dimarray);

#endif
