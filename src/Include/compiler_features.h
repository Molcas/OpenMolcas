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

/* Allocate on assignment */
#if ( __SUNPRO_F90 )
#undef ALLOC_ASSIGN
#else
#define ALLOC_ASSIGN
#endif

/* Pointer bounds remapping */
#if ( __SUNPRO_F90 )
#undef POINTER_REMAP
#else
#define POINTER_REMAP
#endif

/* Trailing zeros in the binary representation (trailz) */
#if ( __PGI )
#undef TRAILING_ZEROS
#else
#define TRAILING_ZEROS
#endif

/* Bit extraction (ibits) with zero length */
#if (( __PGI ) && ( __PGIC__ < 20 ))
#undef IBITS_LEN_ZERO
#else
#define IBITS_LEN_ZERO
#endif

/* c_ptr binding */
#if (( NAGFOR ) && ( __NAG_COMPILER_RELEASE < 61 ) )
#undef C_PTR_BINDING
#else
#define C_PTR_BINDING
#endif

/* Internal procedures as arguments.
With PGI 20 ( __PGIC__ >= 20 ) it compiles, but it appears to be buggy at runtime! */
#if (( __SUNPRO_F90 ) || ( __PGI ))
#undef INTERNAL_PROC_ARG
#else
#define INTERNAL_PROC_ARG
#endif

/* Allows files with no compilable instructions */
#if (( NAGFOR ) || ( __PGI ))
#undef EMPTY_FILES
#else
#define EMPTY_FILES
#endif

/* Storage_size in initialization */
#if (( __GNUC__ ) && ( GCC_VERSION < 70000 ))
#undef SIZE_INITIALIZATION
#else
#define SIZE_INITIALIZATION
#endif

/* Intrinsic functions in initialization */
#if ( __PGI )
#undef INTRINSIC_INITIALIZATION
#else
#define INTRINSIC_INITIALIZATION
#endif

/* Safe character member initialization */
#if (( __GNUC__ ) && ( GCC_VERSION < 80000 ))
#undef CHAR_MEMBER_INIT
#else
#define CHAR_MEMBER_INIT
#endif

/* Empty user-defined type initialization (annoyance in newer ifort) */
#if (( __INTEL_COMPILER ) && ( __INTEL_COMPILER_BUILD_DATE > 20220000 ))
#undef EMPTY_TYPE_INIT
#else
#define EMPTY_TYPE_INIT
#endif
