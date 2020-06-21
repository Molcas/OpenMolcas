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
!***********************************************************************

module constants
    use iso_fortran_env, only: int32, int64, real32, real64
    implicit none
    private
    public :: wp, MPIArg, HDF5Arg
    public :: r4, r8, i4, i8


    ! This is the type of MPI arguments
    ! NOTE: If legacy integer*4 declarations are replaced with integer(MPIArg)
    !       we can support 32bit and 64bit versions.
    !       Which will require appropiate compile flags here.
    integer, parameter :: MPIArg = int32

    ! This is the type of HDF5 arguments
    ! NOTE: If legacy integer*4 declarations are replaced with integer(HDF5Arg)
    !       we can support 32bit and 64bit versions.
    !       Which will require appropiate compile flags here.
    integer, parameter :: HDF5Arg = int32

    ! This is the working precision and should be preferably used.
    integer, parameter :: wp = real64


    ! Although the constants from iso_fortran_env or `selected_real_kind`
    ! are preferred over non-standard real*8 etc.
    ! We define some kinds to refer to the non-standard notation.
    ! **DON'T USE THESE UNLESS YOU EXPLICILTY WANT TO REFER TO real*8 etc.**
    ! `wp` etc. are always preferred.

    real*4 :: r4_example
    real*8 :: r8_example

    integer*4 :: i4_example
    integer*8 :: i8_example

    integer, parameter :: &
        r4 = kind(r4_example), &
        r8 = kind(r8_example), &
        i4 = kind(i4_example), &
        i8 = kind(i8_example)

end module
