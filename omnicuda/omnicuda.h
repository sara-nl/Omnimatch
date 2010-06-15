/**************************************************************************
 * Copyright (C) 2010 Pieter van Beek
 *
 * This file is part of SARAFFT.
 *
 * SARAFFT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SARAFFT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SARAFFT.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

#ifndef OMNICUDA_H
#define OMNICUDA_H

#include <cufft.h>
#include <sarafft.h>

#ifdef __cplusplus
extern "C" {
#endif

void sarafft_init();

sararfftnd_plan sararfft3d_create_plan(
  int nx, int ny, int nz, sarafft_direction dir
);

void sararfftnd_one_real_to_complex(
  sararfftnd_plan p, sarafft_real *data
);

void sararfftnd_one_complex_to_real(
  sararfftnd_plan p, sarafft_complex *data
);

void sararfftnd_destroy_plan(
  sararfftnd_plan plan
);

#ifdef __cplusplus
}
#endif

#endif
