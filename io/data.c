#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "em.h"
#include "mrc.h"

#define EXT_LEN 4

char extension[EXT_LEN]; 

void (*create_output_func)(const char*, int*) = NULL;
void (*read_data_func)(const char*, struct gen_data*) = NULL;
void (*read_header_func)(const char*, struct gen_data*) = NULL;
void (*read_subregion_func)(const char*, struct gen_data*, int*, int*) = NULL;
void (*write_data_func)(const char*, struct gen_data*) = NULL;
void (*write_subregion_func)(const char*, struct gen_data*, int*, int*) = NULL;
void (*write_subsubregion_func)(const char*, struct gen_data*, int*, int*, int*, int*) = NULL;

int findFormat(const char *filename) {
  int retVal = UNKNOWN_FORMAT;
  char *format = 0;
  const char *lastDotPos = strrchr( filename, '.');
  format = strdup(lastDotPos+1);
  if (format != 0) {
    if ( ( strcmp(format, "em") == 0 && isEM( filename ) ) ) {
      strncpy(extension, format, EXT_LEN);
      create_output_func = &create_em;
      read_data_func = &read_em;
      read_header_func = &read_em_header;
      read_subregion_func = &read_em_subregion;
      write_data_func = &write_em;
      write_subregion_func = &write_em_subregion;
      write_subsubregion_func = &write_em_subsubregion;
      retVal = EM_FORMAT;
    } else if ( ( strcmp(format, "mrc") == 0 && isMRC( filename ) ) ) {
      strncpy(extension, format, EXT_LEN);
      create_output_func = &create_mrc;
      read_data_func = &read_mrc;
      read_header_func = &read_mrc_header;
      read_subregion_func = &read_mrc_subregion;
      write_data_func = &write_mrc;
      write_subregion_func = &write_mrc_subregion;
      write_subsubregion_func = &write_mrc_subsubregion;
      retVal = MRC_FORMAT;
    }
    free( format );
  }
  return retVal;
}

void create_output (const char *outfile, int *nr) {
  (*create_output_func)(outfile, nr);
}

void read_data (const char *infile, struct gen_data *data) {
  (*read_data_func)(infile, data);
}

void read_header (const char *infile, struct gen_data *data) {
  (*read_header_func)(infile, data);
}

void read_subregion (const char *infile, struct gen_data *data, int *nr, int *area) {
  (*read_subregion_func)(infile, data, nr, area);
}

void write_data (const char *outfile, struct gen_data *data) {
  (*write_data_func)(outfile, data);
}

void write_subregion (const char *outfile, struct gen_data *data, int *nr, int *area) {
  (*write_subregion_func)(outfile, data, nr, area);
}

void write_subsubregion (const char *outfile, struct gen_data *data, int *nr, int *area,int *offset,int *dimarray) {
  (*write_subsubregion_func)(outfile, data, nr, area, offset, dimarray);
}
