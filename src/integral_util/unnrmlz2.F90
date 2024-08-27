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

subroutine UnNrmlz2(Exp,nPrim,Coeff,nCntrc,iang)

use Constants, only: Two, Three, Four, TwoP34

implicit none
integer nPrim, nCntrc, iAng
real*8 exp(nPrim), Coeff(nPrim,nCntrc)
integer i, j

do i=1,nCntrc
  do j=1,nPrim
    Coeff(j,i) = Coeff(j,i)*(TwoP34*(Four*exp(j))**((Two*iAng+Three)/Four))
  end do
end do

end subroutine UnNrmlz2
