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

module definitions
    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64, input_unit, output_unit
    implicit none
    private
    public :: wp, iwp, MPIInt, HDF5Int
    public :: int32, int64, real32, real64
    public :: i4, i8, r4, r8
    public :: u5, u6

    ! This is the working precision and should be preferably used
    ! (we assume logical kinds are the same as integer kinds).
#ifdef _I8_
    integer(kind=int64), parameter :: iwp = int64
#else
    integer(kind=int32), parameter :: iwp = int32
#endif
    integer(kind=iwp), parameter :: wp = real64

    ! This is the type of MPI arguments
    ! NOTE: If legacy integer*4 declarations are replaced with integer(MPIInt)
    !       we can support 32bit and 64bit versions.
    !       Which will require appropiate compile flags here.
    integer(kind=iwp), parameter :: MPIInt = int32

    ! This is the type of HDF5 arguments
    ! NOTE: If legacy integer*4 declarations are replaced with integer(HDF5Int)
    !       we can support 32bit and 64bit versions.
    !       Which will require appropiate compile flags here.
    integer(kind=iwp), parameter :: HDF5Int = int32

    ! Output and input units, typically 5 & 6, but they could be something else
    integer(kind=iwp), parameter :: u5 = input_unit, &
                                    u6 = output_unit

    ! Although the constants from iso_fortran_env or `selected_real_kind`
    ! are preferred over non-standard real*8 etc.
    ! We define some kinds to refer to the non-standard notation.
    ! **DON'T USE THESE UNLESS YOU EXPLICILTY WANT TO REFER TO real*8 etc.**
    ! `wp` etc. are always preferred.

    real*4 :: r4_example
    real*8 :: r8_example

    integer*4 :: i4_example
    integer*8 :: i8_example

    integer(kind=iwp), parameter :: &
        r4 = kind(r4_example), &
        r8 = kind(r8_example), &
        i4 = kind(i4_example), &
        i8 = kind(i8_example)

end module definitions
