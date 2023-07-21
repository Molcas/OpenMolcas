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

! Stuff for decomposing (ai|bj) integrals or amplitudes in MP2:
module chomp2_dec

public :: InCore, iOption_MP2CD, NowSym
logical InCore(8)
integer iOption_MP2CD, NowSym
real*8, pointer, contiguous :: EOcc(:) => null(), EVir(:) => null()

end module chomp2_dec
