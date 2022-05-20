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
! Object: to setup the roots and weights of the Hermite polynomials    *
!         for the evaluation of one electron integrals.                *
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
integer(kind=iwp) :: nDiff
integer(kind=iwp) :: i_0000, i_1111, i_2222, i_3333, IDEG, IDH, iHer, IR, IROOT, IW, j, j_0000, j_1111, j_2222, j_3333, K, n_1111, &
                     n_2222, nMem
real(kind=wp) :: Alpha, B, b_1111, c_0000, CORR, DELTA, HDER, R, RSUM, W, w_1111, w_2222, w_3333, w_4444, w_5555, X, Z
real(kind=wp), allocatable :: Beta(:), BInv(:), Herm(:)

! 1) Hermite-Gauss
! 2) Rys-Gauss (asymtotic formula)

n_1111 = (2*S%iAngMx+nPrp+2+nDiff)/2
n_2222 = 4*S%iAngMx+2+nDiff

if (allocated(HerR) .and. (max(n_1111,n_2222) <= MaxHer)) then
  return
else if (allocated(HerR)) then
  call Free_HerRW()
end if
MaxHer = max(n_1111,n_2222)
call mma_allocate(iHerR,MaxHer,label='iHerR')
call mma_allocate(iHerW,MaxHer,label='iHerW')

! Set up square of roots and weights for Hermite polynomials

nMem = (MaxHer*MaxHer+MaxHer)/2
call mma_Allocate(HerR,nMem,label='HerR')
iHerR(1) = 1
call dCopy_(nMem,[Zero],0,HerR,1)
call mma_allocate(HerW,nMem,label='HerW')
iHerW(1) = 1
call dCopy_(nMem,[Zero],0,HerW,1)
call mma_allocate(Beta,MaxHer,label='Beta')
call dCopy_(MaxHer,[Zero],0,Beta,1)
call mma_allocate(BInv,MaxHer,label='BInv')
call dCopy_(MaxHer,[Zero],0,BInv,1)
call mma_allocate(Herm,MaxHer+1,label='Herm')
call dCopy_(MaxHer+1,[Zero],0,Herm,1)
do K=1,MaxHer
  b_1111 = HALF*real(K,kind=wp)
  B = sqrt(b_1111)
  Beta(K) = B
  BInv(K) = One/B
end do
HerR(iHerR(1)) = Zero
HerR(iHerR(1)+2) = sqrt(HALF)
HerR(iHerR(1)+1) = -HerR(iHerR(1)+2)
HerW(iHerW(1)) = sqrt(PI)
HerW(iHerW(1)+1) = HerW(iHerW(1))*HALF
HerW(iHerW(1)+2) = HerW(iHerW(1)+1)
Herm(1) = One/sqrt(HerW(iHerW(1)))
do iHer=2,MaxHer
  i_1111 = (iHer*iHer-iHer)/2
  iHerR(iHer) = iHerR(1)+i_1111
  iHerW(iHer) = iHerW(1)+i_1111
end do

Alpha = BInv(1)
do IDEG=3,MaxHer
  i_0000 = (IDEG*IDEG-IDEG)/2
  IR = iHerR(1)-1+i_0000
  IW = iHerW(1)-1+i_0000
  IDH = IDEG/2
  i_1111 = IR+IDH+1
  i_3333 = i_1111-IDEG
  w_3333 = HerR(i_3333)
  i_2222 = i_3333+1
  w_2222 = HerR(i_2222)
  X = (w_2222-w_3333)/Two
  HerR(i_1111) = Zero
  do IROOT=2,IDEG,2
    j_0000 = IROOT/2
    j_1111 = IR+j_0000
    j_2222 = IR-IDEG+1+j_0000
    j_3333 = IR+IDEG+1-j_0000
    R = HerR(j_2222)-X
    HerR(j_1111) = R
    HerR(j_3333) = -R
  end do
  do IROOT=1,IDH
    j_0000 = IR+IROOT
    Z = HerR(j_0000)
    CORR = Zero
    do j=1,ideg
      if (j /= iroot) then
        c_0000 = Z-HerR(IR+J)
        CORR = CORR+(One/c_0000)
      end if
    end do
    do
      Herm(2) = Z*Herm(1)*Alpha
      do K=1,IDEG-1
        w_1111 = Herm(K+1)
        w_3333 = Herm(K)
        w_4444 = Beta(K)
        w_5555 = BInv(K+1)
        w_2222 = (Z*w_1111-w_4444*w_3333)*w_5555
        Herm(K+2) = w_2222
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
    j_0000 = IR+IROOT
    Z = HerR(j_0000)
    Herm(2) = Z*Herm(1)*Alpha
    RSUM = Herm(1)**2
    RSUM = RSUM+Herm(2)**2
    do K=1,IDEG-2
      w_1111 = Herm(K+1)
      w_3333 = Herm(K)
      w_4444 = Beta(K)
      w_5555 = BInv(K+1)
      w_2222 = (Z*w_1111-w_4444*w_3333)*w_5555
      Herm(K+2) = w_2222
      RSUM = RSUM+w_2222*w_2222
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
