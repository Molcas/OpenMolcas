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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
***********************************************************************/

/*
This preprocessor include file may help finding out if some features
are available in the current compiler and enable workarounds if needed.

Minimal versions are probably not accurate, and compiler support is
incomplete.
*/

#if (__GNUC__)
#define GCC_VERSION (__GNUC__ * 10000  + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#else
#define GCC_VERSION 0
#endif

/* Allocatable characters */
#if (__GNUC__ && GCC_VERSION < 40800 )
#undef ALLOC_CHAR
#else
#define ALLOC_CHAR
#endif

/* Allocate on assignment */
#if ((!defined(ALLOC_CHAR)) || \
     ( __SUNPRO_F90 ))
#undef ALLOC_ASSIGN
#else
#define ALLOC_ASSIGN
#endif

/* Allocatable scalars */
#if (__GNUC__ && GCC_VERSION < 40500 )
#undef ALLOC_SCAL
#else
#define ALLOC_SCAL
#endif

/* Pointer bounds remapping */
#if (__GNUC__ && GCC_VERSION < 40700 )
#undef POINTER_REMAP
#else
#define POINTER_REMAP
#endif

/* Parity of the binary representation (poppar) */
#if ((__GNUC__ && GCC_VERSION < 40600 ) || \
     ( __INTEL_COMPILER && __INTEL_COMPILER < 1300 ))
#undef BINARY_PARITY
#else
#define BINARY_PARITY
#endif

/* Trailing zeros in the binary representation (trailz) */
#if ((__GNUC__ && GCC_VERSION < 40600 ) || \
     ( __INTEL_COMPILER && __INTEL_COMPILER < 1300 ) || \
     ( __PGI ))
#undef TRAILING_ZEROS
#else
#define TRAILING_ZEROS
#endif

/* c_ptr binding */
#if (NAGFOR && __NAG_COMPILER_RELEASE < 61 )
#undef C_PTR_BINDING
#else
#define C_PTR_BINDING
#endif
