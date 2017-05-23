/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Per-Olof Widmark                                       *
***********************************************************************/
/**************************************************************************/
/*                                                                        */
/* This routine copy a number of bytes.                                   */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/*                                                                        */
/**************************************************************************/
#include "molcastype.h"
#ifndef _WIN32
#include <strings.h>
#else
#include <string.h>
#endif
#ifdef _CAPITALS_
#define bytecopy  BYTECOPY
#else
#ifndef ADD_
#define bytecopy  bytecopy_
#endif
#endif
#ifdef __cplusplus
extern "C" {
#endif
void bytecopy(unsigned char *from, unsigned char *to, INT *nbytes) {
#ifndef _WIN32
     bcopy(from,to,(size_t)(*nbytes));
#else
     memcpy(to,from,(size_t)(*nbytes));
#endif
}
#ifdef __cplusplus
}
#endif
