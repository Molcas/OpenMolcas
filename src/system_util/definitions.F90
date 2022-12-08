!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2020, Oskar Weser                                      *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************

module Definitions

use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64, error_unit, input_unit, output_unit
use, intrinsic :: iso_c_binding, only: c_double, c_int, c_size_t
#   ifdef _I8_
use, intrinsic :: iso_c_binding, only: c_long
#   endif

implicit none
private

public :: wp, iwp, byte, DefInt, MPIInt, HDF5Int
public :: BLASInt, BLASR4, BLASR8, CUDAInt
public :: LibxcInt, LibxcReal, LibxcSize
public :: MOLCAS_C_INT, MOLCAS_C_REAL
public :: i1, i4, i8, r4, r8
public :: ItoB, RtoB, RtoI, CtoR
public :: u0, u5, u6

! This is the working precision and should be preferably used
! (we assume logical kinds are the same as integer kinds).
#ifdef _I8_
integer(kind=int64), parameter :: iwp = int64, MOLCAS_C_INT = c_long
#else
integer(kind=int32), parameter :: iwp = int32, MOLCAS_C_INT = c_int
#endif
integer(kind=iwp), parameter :: wp = real64, MOLCAS_C_REAL = c_double

! Single byte integer
integer(kind=iwp), parameter :: byte = int8

! "default" integer, without using `-i8` flag or equivalent,
! this is needed for some intrinsic calls in some compilers
integer(kind=iwp), parameter :: DefInt = int32

! Types for BLAS/LAPACK calls
#if defined(LINALG_I4) && defined(_I8_)
integer(kind=iwp), parameter :: BLASInt = int32
#else
integer(kind=iwp), parameter :: BLASInt = iwp
#endif
integer(kind=iwp), parameter :: BLASR4 = real32, BLASR8 = real64, CUDAInt = int32

! This is the type of Libxc arguments
integer(kind=iwp), parameter :: LibxcInt = c_int, &
                                LibxcReal = c_double, &
                                LibxcSize = c_size_t

! This is the type of MPI arguments
! NOTE: If legacy `integer*4` declarations are replaced with integer(MPIInt)
!       we can support 32bit and 64bit versions.
!       Which will require appropiate compile flags here.
integer(kind=iwp), parameter :: MPIInt = int32

! This is the type of HDF5 arguments
! NOTE: If legacy `integer*4` declarations are replaced with integer(HDF5Int)
!       we can support 32bit and 64bit versions.
!       Which will require appropiate compile flags here.
integer(kind=iwp), parameter :: HDF5Int = int32

! Size ratios between the different types
! (assume 'a' uses 1 byte)
#include "compiler_features.h"
integer(kind=iwp), parameter :: &
#ifdef SIZE_INITIALIZATION
                                ItoB = storage_size(1_iwp)/storage_size('a'), &
                                RtoB = storage_size(1.0_wp)/storage_size('a'), &
                                RtoI = storage_size(1.0_wp)/storage_size(1_iwp), &
                                CtoR = storage_size((1.0_wp,0.0_wp))/storage_size(1.0_wp)
#elif defined(_I8_)
                                ItoB = 8, RtoB = 8, RtoI = 1, CtoR = 2
#else
                                ItoB = 4, RtoB = 8, RtoI = 2, CtoR = 2
#endif

! Output, input and error units, typically 6, 5 & 0, but they could be something else
integer(kind=iwp), parameter :: u0 = error_unit, &
                                u5 = input_unit, &
                                u6 = output_unit

! Although the constants from `iso_fortran_env` or `selected_real_kind`
! are preferred over non-standard `real*8` etc.
! We define some kinds to refer to the non-standard notation.
! **DON'T USE THESE UNLESS YOU EXPLICILTY WANT TO REFER TO `real*8` etc.**
! `wp` etc. are always preferred.

real*4 :: r4_example
real*8 :: r8_example

integer*1 :: i1_example
integer*4 :: i4_example
integer*8 :: i8_example

integer(kind=iwp), parameter :: r4 = kind(r4_example), &
                                r8 = kind(r8_example), &
                                i1 = kind(i1_example), &
                                i4 = kind(i4_example), &
                                i8 = kind(i8_example)

end module Definitions
