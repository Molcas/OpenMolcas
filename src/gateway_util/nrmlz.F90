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

subroutine Nrmlz(rExp,nPrim,Coeff,nCntrc,iAng)
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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, OneHalf, TwoP34
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPrim, nCntrc, iAng
real(kind=wp), intent(in) :: rExp(nPrim)
real(kind=wp), intent(inout) :: Coeff(nPrim,nCntrc)
integer(kind=iwp) :: i, iExp, j, jExp, nScrt1, nScrt2
real(kind=wp) :: Power, Pro_ij, Qtemp, Rtemp, Sum_ij, Temp, vR2, vRR
real(kind=wp), allocatable :: Scrt1(:), Scrt2(:)

if (nPrim*nCntrc == 0) return

nScrt1 = nPrim**2
nScrt2 = nPrim*nCntrc
call mma_allocate(Scrt1,nScrt1)
call mma_allocate(Scrt2,nScrt2)

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
    Pro_ij = sqrt(rExp(iExp)*rExp(jExp))
    Sum_ij = (rExp(iExp)+rExp(jExp))/Two
    Power = real(iAng,kind=wp)+OneHalf
    Temp = (Pro_ij/Sum_ij)**Power
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
#ifdef _DEBUGPRINT_
call RecPrt(' Overlap Contracted',' ',Scrt1,nCntrc,1)
#endif

! Normalize coefficients, i.e. combine the normalization factor
! of the primitive and the overlap of the unnormalized contracted
! basis function.

do i=1,nCntrc
  if (abs(Scrt1(i)) < 1.0e-12_wp) then
    call WarningMessage(2,'; Error in contraction matrix, zero column; ; Abend in subroutine NRMLZ')
    call Abend()
  end if
end do

Rtemp = Half*real(iAng,kind=wp)+0.75_wp
Qtemp = Two**(iAng+1)*sqrt(Two)*TwoP34
do i=1,nCntrc
  vRR = Scrt1(i)**(-Half)
  do j=1,nPrim
    vR2 = rExp(j)**Rtemp
    Coeff(j,i) = Coeff(j,i)*Qtemp*vRR*vR2
  end do
end do
if ((nPrim == 1) .and. (nCntrc == 1) .and. (rExp(1) == Zero)) then
  Coeff(1,1) = One
end if
#ifdef _DEBUGPRINT_
call Recprt(' In Nrmlz: Normalized coefficients',' ',Coeff,nPrim,nCntrc)
#endif

call mma_deallocate(Scrt2)
call mma_deallocate(Scrt1)

return

end subroutine Nrmlz
