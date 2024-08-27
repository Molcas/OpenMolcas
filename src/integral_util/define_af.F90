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

private

#ifndef _DEMO_
integer, parameter :: iTabMx = 15
character(len=1), parameter :: AngTp(0:iTabMx) = ['s','p','d','f','g','h','i','k','l','m','n','o','q','r','t','u']
#else
integer, parameter :: iTabMx = 3
character(len=1), parameter :: AngTp(0:iTabMx) = ['s','p','d']
#endif

real*8 Binom(0:2*iTabMx,-1:2*iTabMx)
integer, parameter :: lab = 2*iTabMx+1, ipMax = lab*(lab+1)*(lab+2)/6
integer iCan(3,ipMax)

public :: iTabMx, AngTp, Binom, lab, ipMax, iCan

end module Define_af
