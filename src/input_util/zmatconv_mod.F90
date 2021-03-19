!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module ZMatConv_Mod

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MaxAtoms = 1000
character(len=48), allocatable :: Base(:)      ! Base: Basis Set string
logical(kind=iwp), allocatable :: BasAva(:), & ! Atom with available Basis Set
                                  BasReq(:)    ! Atom requiring Basis Set
integer(kind=iwp), allocatable :: NAT(:), &    ! Atomic number
                                  iZmat(:,:)   ! Z-Mat indices
real(kind=wp), allocatable :: Zmat(:,:), &     ! Z-Mat coordinates
                              Coords(:,:)      ! Atomic coordinates
character(len=5), allocatable :: Symbols(:)    ! Atomic symbol with index. (C12 or Hn)

public :: Base, BasAva, BasReq, Symbols, MaxAtoms, NAT, iZmat, Zmat, Coords

end module ZMatConv_Mod
