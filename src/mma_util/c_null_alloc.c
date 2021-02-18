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
This function can be used to mark a Fortran allocatable array as deallocated.
Note that it doesn't deallocate any memory, it just overwrites the array
descriptor. It is useful for resetting allocatable components in a derived
type array, after the contents have been garbled.
*/
#include <stdint.h>
#ifdef _CAPITALS_
#define c_null_alloc C_NULL_ALLOC
#else
#ifndef ADD_
#define c_null_alloc c_null_alloc_
#endif
#endif
void c_null_alloc(intptr_t* A) { *A = 0; }
