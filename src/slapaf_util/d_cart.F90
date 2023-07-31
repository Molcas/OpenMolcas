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

function D_Cart(Ind,nSym)

use Slapaf_Info, only: nStab
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: D_Cart
integer(kind=iwp), intent(in) :: Ind, nSym
integer(kind=iwp) :: iAtom, iDeg, nU_A

!                                                                      *
!***********************************************************************
!                                                                      *
! Cartesian coordinate  (iAtom)

D_Cart = Zero

iAtom = Ind

nU_A = nStab(iAtom)

! Now evaluate the degeneracy of the cartesian

iDeg = nSym/nU_A
d_cart = real(iDeg,kind=wp)

!write(u6,*) ' d_cart=',d_cart

return

end function D_Cart
