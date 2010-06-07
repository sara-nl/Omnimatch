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

#ifndef SARAFFT_H
#define SARAFFT_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef USE_GPUS

#else // #ifdef USE_GPUS

#endif // #ifdef USE_GPUS

#include <sfftw.h>
#include <srfftw.h>

typedef fftw_real sarafft_real;
typedef fftw_complex sarafft_complex;

#ifdef __cplusplus
} // extern "C"
#endif

#endif // #ifndef SARAFFT_H
