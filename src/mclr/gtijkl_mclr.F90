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

function GTIJKL_MCLR(I,J,K,L)
! Obtain  integral (I J ! K L)
! where I,J,K and l refers to active orbitals in type ordering

use Index_Functions, only: iTri
use MCLR_Data, only: Int2, IREOTS
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GTIJKL_MCLR
integer(kind=iwp), intent(in) :: I, J, K, L
integer(kind=iwp) :: i_Abs, ij, j_Abs, k_Abs, kl, l_Abs

I_ABS = IREOTS(I)

J_ABS = IREOTS(J)

K_ABS = IREOTS(K)

L_ABS = IREOTS(L)

IJ = iTri(I_ABS,J_ABS)
KL = iTri(K_ABS,L_ABS)

GTIJKL_MCLR = INT2(iTri(IJ,KL))

end function GTIJKL_MCLR
