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

real*8 function GTIJKL_MCLR(I,J,K,L)
! Obtain  integral (I J ! K L)
! where I,J,K and l refers to active orbitals in type ordering

use Index_Functions, only: iTri
use MCLR_Data, only: Int2
use MCLR_Data, only: IREOTS

implicit none
integer, intent(in) :: I, J, K, L
integer iAbs, jAbs, kAbs, lAbs, ij, kl

IABS = IREOTS(I)

JABS = IREOTS(J)

KABS = IREOTS(K)

LABS = IREOTS(L)

IJ = iTri(IABS,JABS)
KL = iTri(KABS,LABS)

GTIJKL_MCLR = INT2(iTri(IJ,KL))

end function GTIJKL_MCLR
