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

use Her_RW
use Sizes_of_Seward, only: S

implicit real*8(A-H,O-Z)
#include "stdalloc.fh"
#include "real.fh"
#include "status.fh"
real*8, dimension(:), allocatable :: Beta, BInv, Herm

if (nPrp > nPrpMx) then
  write(6,*) 'nPrp, nPrpMx=',nPrp,nPrpMx
  call WarningMessage(2,'SetHer: nPrp too large!')
  call Abend()
end if

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
call dCopy_(nMem,[0.0d0],0,HerR,1)
call mma_allocate(HerW,nMem,label='HerW')
iHerW(1) = 1
call dCopy_(nMem,[0.0d0],0,HerW,1)
call mma_allocate(Beta,MaxHer,label='Beta')
call dCopy_(MaxHer,[0.0d0],0,Beta,1)
call mma_allocate(BInv,MaxHer,label='BInv')
call dCopy_(MaxHer,[0.0d0],0,BInv,1)
call mma_allocate(Herm,MaxHer+1,label='Herm')
call dCopy_(MaxHer+1,[0.0d0],0,Herm,1)
do K=1,MaxHer
  b_1111 = HALF*dble(K)
  B = sqrt(b_1111)
  Beta(K) = B
  BInv(K) = 1.0d0/B
end do
HerR(iHerR(1)) = 0.0d0
HerR(iHerR(1)+2) = sqrt(HALF)
HerR(iHerR(1)+1) = -HerR(iHerR(1)+2)
HerW(iHerW(1)) = sqrt(PI)
HerW(iHerW(1)+1) = HerW(iHerW(1))*HALF
HerW(iHerW(1)+2) = HerW(iHerW(1)+1)
Herm(1) = 1.0d0/sqrt(HerW(iHerW(1)))
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
  X = (w_2222-w_3333)/2.0d0
  HerR(i_1111) = 0.0d0
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
    CORR = 0.0d0
    do j=1,ideg
      if (j /= iroot) then
        c_0000 = Z-HerR(IR+J)
        CORR = CORR+(1.0d0/c_0000)
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
      HDER = 2.0d0*Beta(IDEG)*Herm(IDEG)
      DELTA = -Herm(IDEG+1)/(HDER-CORR*Herm(IDEG+1))
      Z = Z+DELTA
      if (abs(DELTA) <= 1.0d-8) exit
      if (abs(DELTA) > 1.0d8) then
        call WarningMessage(1,'Warning: large value in sether')
        !write(6,*) delta
      end if
    end do
    HerR(IR+IROOT) = Z
    HerR(IR+IDEG+1-IROOT) = -Z
  end do
  do IROOT=1,IDH+1
    j_0000 = IR+IROOT
    Z = HerR(j_0000)
    Herm(2) = Z*Herm(1)*Alpha
    SUM = Herm(1)**2
    SUM = SUM+Herm(2)**2
    do K=1,IDEG-2
      w_1111 = Herm(K+1)
      w_3333 = Herm(K)
      w_4444 = Beta(K)
      w_5555 = BInv(K+1)
      w_2222 = (Z*w_1111-w_4444*w_3333)*w_5555
      Herm(K+2) = w_2222
      SUM = SUM+w_2222*w_2222
    end do
    W = 1.0d0/SUM
    HerW(IW+IROOT) = W
    HerW(IW+IDEG+1-IROOT) = W
  end do
end do
call mma_deallocate(Beta)
call mma_deallocate(BInv)
call mma_deallocate(Herm)

!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call TriPrt(' Hermite roots',' ',HerR(iHerR(1)),MaxHer)
call TriPrt(' Hermite weights',' ',HerW(iHerW(1)),MaxHer)
write(6,*) ' MaxHer=',MaxHer,nPrp,S%iAngMx
#endif

return

end subroutine SetHer
