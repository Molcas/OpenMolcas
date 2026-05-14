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
! This defines the highest angular momentum quantum number which
! Seward will be able to treat.

module Define_af

use Definitions, only: wp, iwp

implicit none
private

#ifndef _DEMO_
integer(kind=iwp), parameter :: iTabMx = 15
character, parameter :: AngTp(0:iTabMx) = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']
#else
integer(kind=iwp), parameter :: iTabMx = 3
character, parameter :: AngTp(0:iTabMx) = ['s','p','d']
#endif

integer(kind=iwp), parameter :: lab = 2*iTabMx+1, ipMax = lab*(lab+1)*(lab+2)/6
real(kind=wp) :: Binom(0:2*iTabMx,-1:2*iTabMx)
integer(kind=iwp) :: iCan(3,ipMax)

public :: AngTp, Binom, iCan, ipMax, iTabMx, lab

end module Define_af
