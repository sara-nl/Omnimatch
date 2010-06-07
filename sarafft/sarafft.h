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

#ifdef __cplusplus
}
#endif

#endif // #ifndef SARAFFT_H