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
! Copyright (C) 1992, Per Ake Malmqvist                                *
!               1992, Roland Lindh                                     *
!***********************************************************************

subroutine SetHer(nDiff)
!***********************************************************************
!                                                                      *
! Object: to set up the roots and weights of the Hermite polynomials   *
!         for the evaluation of one-electron integrals.                *
!                                                                      *
!    Authors: Per-AAke Malmqvist and Roland Lindh,                     *
!             March 1992.                                              *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW, MaxHer, nPrp
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Pi
use Definitions, only: wp, iwp
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nDiff
integer(kind=iwp) :: i_0, i_1, IDEG, IDH, iHer, IR, IROOT, IW, j, j_0, K, n_1, n_2, nMem
real(kind=wp) :: Alpha, CORR, DELTA, HDER, R, RSUM, W, X, Z
real(kind=wp), allocatable :: Beta(:), BInv(:), Herm(:)

! 1) Hermite-Gauss
! 2) Rys-Gauss (asymptotic formula)

n_1 = (2*S%iAngMx+nPrp+2+nDiff)/2
n_2 = 4*S%iAngMx+4+nDiff

if (allocated(HerR) .and. (max(n_1,n_2) <= MaxHer)) then
  return
else if (allocated(HerR)) then
  call Free_HerRW()
end if
MaxHer = max(n_1,n_2)
call mma_allocate(iHerR,MaxHer,label='iHerR')
iHerR(1) = 1
call mma_allocate(iHerW,MaxHer,label='iHerW')
iHerW(1) = 1

! Set up square of roots and weights for Hermite polynomials

nMem = (MaxHer*MaxHer+MaxHer)/2
call mma_Allocate(HerR,nMem,label='HerR')
HerR(:) = Zero
call mma_allocate(HerW,nMem,label='HerW')
HerW(:) = Zero
call mma_allocate(Beta,MaxHer,label='Beta')
call mma_allocate(BInv,MaxHer,label='BInv')
call mma_allocate(Herm,MaxHer+1,label='Herm')
Herm(:) = Zero
do K=1,MaxHer
  Beta(K) = sqrt(Half*K)
end do
BInv(:) = One/Beta
HerR(iHerR(1)) = Zero
HerR(iHerR(1)+2) = sqrt(HALF)
HerR(iHerR(1)+1) = -HerR(iHerR(1)+2)
HerW(iHerW(1)) = sqrt(PI)
HerW(iHerW(1)+1) = HerW(iHerW(1))*HALF
HerW(iHerW(1)+2) = HerW(iHerW(1)+1)
Herm(1) = One/sqrt(HerW(iHerW(1)))
do iHer=2,MaxHer
  i_1 = (iHer*iHer-iHer)/2
  iHerR(iHer) = iHerR(1)+i_1
  iHerW(iHer) = iHerW(1)+i_1
end do

Alpha = BInv(1)
do IDEG=3,MaxHer
  i_0 = (IDEG*IDEG-IDEG)/2
  IR = iHerR(1)-1+i_0
  IW = iHerW(1)-1+i_0
  IDH = IDEG/2
  i_1 = IR+IDH+1
  X = Half*(HerR(i_1-IDEG+1)-HerR(i_1-IDEG))
  HerR(i_1) = Zero
  do IROOT=2,IDEG,2
    j_0 = IROOT/2
    R = HerR(IR-IDEG+1+j_0)-X
    HerR(IR+j_0) = R
    HerR(IR+IDEG+1-j_0) = -R
  end do
  do IROOT=1,IDH
    Z = HerR(IR+IROOT)
    CORR = Zero
    do j=1,ideg
      if (j /= iroot) CORR = CORR+(One/Z-HerR(IR+J))
    end do
    do
      Herm(2) = Z*Herm(1)*Alpha
      do K=1,IDEG-1
        Herm(K+2) = (Z*Herm(K+1)-Beta(K)*Herm(K))*BInv(K+1)
      end do
      HDER = Two*Beta(IDEG)*Herm(IDEG)
      DELTA = -Herm(IDEG+1)/(HDER-CORR*Herm(IDEG+1))
      Z = Z+DELTA
      if (abs(DELTA) <= 1.0e-8_wp) exit
      if (abs(DELTA) > 1.0e8_wp) then
        call WarningMessage(1,'Warning: large value in sether')
        !write(u6,*) delta
      end if
    end do
    HerR(IR+IROOT) = Z
    HerR(IR+IDEG+1-IROOT) = -Z
  end do
  do IROOT=1,IDH+1
    Z = HerR(IR+IROOT)
    Herm(2) = Z*Herm(1)*Alpha
    RSUM = Herm(1)**2
    RSUM = RSUM+Herm(2)**2
    do K=1,IDEG-2
      Herm(K+2) = (Z*Herm(K+1)-Beta(K)*Herm(K))*BInv(K+1)
      RSUM = RSUM+Herm(K+2)**2
    end do
    W = One/RSUM
    HerW(IW+IROOT) = W
    HerW(IW+IDEG+1-IROOT) = W
  end do
end do
call mma_deallocate(Beta)
call mma_deallocate(BInv)
call mma_deallocate(Herm)

#ifdef _DEBUGPRINT_
call TriPrt(' Hermite roots',' ',HerR(iHerR(1)),MaxHer)
call TriPrt(' Hermite weights',' ',HerW(iHerW(1)),MaxHer)
write(u6,*) ' MaxHer=',MaxHer,nPrp,S%iAngMx
#endif

return

end subroutine SetHer
