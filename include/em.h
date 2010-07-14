/*
 *  em.h : header for EM file read/write functions
 *
 */

#ifndef EM_H
#define EM_H

#include "datatypes.h"

int  isEM( const char *filename);
void create_em (const char *outfile, int *nr);
void read_em (const char *infile, struct gen_data *inemdata);
void read_em_header (const char *infile, struct gen_data *inemdata);
void read_em_subregion (const char *infile, struct gen_data *inemdata, int *nr, int *area);
void write_em (const char *outfile, struct gen_data *outemdata);
void write_em_subregion (const char *infile, struct gen_data *inemdata, int *nr, int *area);
void write_em_subsubregion (const char *outfile, struct gen_data *outemdata, int *nr, int *area,int *offset,int *dimarray);

#endif
