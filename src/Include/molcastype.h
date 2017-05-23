/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
/*
 *
 * header file for int definition
 *
 */

#ifndef _MOLCAS_TYPE_H_
#define _MOLCAS_TYPE_H_

#ifdef _I8_
#ifdef _x86_I8_           /* ---- x86 I8  ---- */
#define INT long long int
#define UN_INT unsigned long long int
#define INT_FORMAT "%lld"
#else
#define INT long int      /* ----   I8   ----- */
#define UN_INT unsigned long int
#define INT_FORMAT "%ld"
#endif
#else
#define INT int           /* ----   I4   ----- */
#define UN_INT unsigned int
#define INT_FORMAT "%d"
#endif
#define LIFMT(X) (long int) (X)

#endif
