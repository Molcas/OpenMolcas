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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine SECORD(A,B,C,FAC,NAL,NBL,NSIJ,IFT)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: A(*), B(*), FAC
real(kind=wp), intent(_OUT_) :: C(*)
integer(kind=iwp), intent(in) :: NAL, NBL, NSIJ, IFT
integer(kind=iwp) :: IAB, NA, NA1, NAA, NB, NBB

IAB = 0
NAA = 0
do NA=1,NAL
  NBB = 0
  NA1 = NBL
  if (NSIJ == 1) NA1 = NA-1
  do NB=1,NA1
    IAB = IAB+1
    if (IFT == 0) C(IAB) = B(NAA+NB)+A(NBB+NA)
    if (IFT == 1) C(IAB) = B(NAA+NB)-A(NBB+NA)
    NBB = NBB+NAL
  end do
  if (NSIJ == 1) then
    IAB = IAB+1
    C(IAB) = FAC*A(NAA+NA)
  end if
  NAA = NAA+NBL
end do

return

end subroutine SECORD
