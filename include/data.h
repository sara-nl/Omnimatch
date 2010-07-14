#ifndef DATA_H
#define DATA_H

#include "datatypes.h"

enum {UNKNOWN_FORMAT=0, EM_FORMAT, MRC_FORMAT};

int findFormat(const char *filename);
void create_output (const char *outfile, int *nr);
void read_data (const char *infile, struct gen_data *data);
void read_header (const char *infile, struct gen_data *data);
void read_subregion (const char *infile, struct gen_data *data, int *nr, int *area);
void write_data (const char *outfile, struct gen_data *data);
void write_subregion (const char *infile, struct gen_data *data, int *nr, int *area);
void write_subsubregion (const char *outfile, struct gen_data *data, int *nr, int *area,int *offset,int *dimarray);

#endif
