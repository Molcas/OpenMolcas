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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine Nrmlx(rExp,nPrim,Coeff,nCntrc,Scrt1,nScrt1,Scrt2,nScrt2,iAng)
!***********************************************************************
!                                                                      *
! Object: normalize the contraction coefficients with respect to the   *
!         radial overlap.                                              *
!                                                                      *
! Called from: Input                                                   *
!                                                                      *
! Calling    : RecPrt                                                  *
!              DGEMM_  (ESSL)                                          *
!              DnDot   (ESSL)                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!***********************************************************************

use Constants, only: Zero, One, Two, Three, Three
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nPrim, nCntrc, nScrt1, nScrt2, iAng
real(kind=wp), intent(in) :: rExp(nPrim)
real(kind=wp), intent(inout) :: Coeff(nPrim,nCntrc)
real(kind=wp), intent(out) :: Scrt1(nScrt1), Scrt2(nScrt2)
integer(kind=iwp) :: i, iExp, jExp
real(kind=wp) :: Temp

#ifdef _DEBUGPRINT_
write(u6,*) ' In Nrmlz: iAng=',iAng
call RecPrt(' In Nrmlz: Coefficients',' ',Coeff,nPrim,nCntrc)
call RecPrt(' In Nrmlz: Exponents',' ',rExp,nPrim,1)
#endif

! Normalize the coefficients (only radial normalization)

! Compute overlap for all primitives of this shell.
! This formula includes the normalization constant for each
! primitive as well as the overlap factor of between the primitives.
! Hence, the overlap matrix elements correspond to those of the
! normalized primitive gaussian functions.

do iExp=1,nPrim
  do jExp=1,iExp-1
    Temp = (Two*sqrt(rExp(iExp)*rExp(jExp))/(rExp(iExp)+rExp(jExp)))**(real(iAng,kind=wp)+Three/Two)
    Scrt1(nPrim*(iExp-1)+jExp) = Temp
    Scrt1(nPrim*(jExp-1)+iExp) = Temp
  end do
  Scrt1(nPrim*(iExp-1)+iExp) = One
end do
! Contract right side
call DGEMM_('N','N',nPrim,nCntrc,nPrim,One,Scrt1,nPrim,Coeff,nPrim,Zero,Scrt2,nPrim)
#ifdef _DEBUGPRINT_
call RecPrt(' Overlap primitives',' ',Scrt1,nPrim,nPrim)
call RecPrt(' Overlap PrimCon',' ',Scrt2,nPrim,nCntrc)
#endif
! Compute the overlap for each contracted basis function, <i|i>
call DnDot(nCntrc,nPrim,Scrt1,1,1,Scrt2,1,nPrim,Coeff,1,nPrim)
! Normalize coefficients, i.e. combine the normalization factor
! of the primitive and the overlap of the unnormalized contracted
! basis function.
do i=1,nCntrc
  Coeff(:,i) = Coeff(:,i)/sqrt(Scrt1(i))
end do
#ifdef _DEBUGPRINT_
call Recprt(' In Nrmlz: Normalized coefficients',' ',Coeff,nPrim,nCntrc)
#endif

return

end subroutine Nrmlx
