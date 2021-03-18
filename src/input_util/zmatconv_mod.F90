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

integer(kind=iwp), parameter :: MaxAtoms = 256
#include "periodic_table.fh"
character(len=48) :: Base(Num_Elem)      ! Base: Basis Set string
logical(kind=iwp) :: BasAva(Num_Elem), & ! Atom with available Basis Set
                     BasReq(Num_Elem)    ! Atom requiring Basis Set
character(len=5) :: Symbols(MaxAtoms)    ! Atomic symbol with index. (C12 or Hn)
integer(kind=iwp) :: NAT(MaxAtoms), &    ! Atomic number
                     iZmat(MaxAtoms,3)   ! Z-Mat indices
real(kind=wp) :: ZMat(MaxAtoms,3), &     ! Z-Mat coordinates
                 Coords(MaxAtoms,3)      ! Atomic coordinates

public :: Num_Elem, PTab, MaxAtoms, Base, BasAva, BasReq, Symbols, NAT, iZmat, ZMat, Coords

end module ZMatConv_Mod
