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

#include "compiler_features.h"
#ifndef _ALLOC_BIND_C_
#ifdef _CAPITALS_
# define c_null_alloc C_NULL_ALLOC
# define c_null_alloc2 C_NULL_ALLOC2
# define c_null_alloc3 C_NULL_ALLOC3
#else
# ifndef ADD_
#   define c_null_alloc c_null_alloc_
#   define c_null_alloc2 c_null_alloc2_
#   define c_null_alloc3 c_null_alloc3_
# endif
#endif
#endif

#include <stdint.h>

void c_null_alloc(intptr_t *A) {
  *A = 0;
}

/* silly wrappers to work around compilers complaining about mismatched interfaces */
void c_null_alloc2(intptr_t *A) {
  c_null_alloc(A);
}

void c_null_alloc3(intptr_t *A) {
  c_null_alloc(A);
}
