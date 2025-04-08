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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)

use ipPage, only: ipget, W
use MCLR_Data, only: nConf1, nAcPr2, ipCI, ipMat, nDens
use MCLR_Data, only: XISPSM
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: nRoots, State_Sym, nSym, nBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp

implicit none
! Output
real*8, dimension(nConf1,nRoots) :: CSFOK
real*8, dimension(nRoots,nRoots) :: LOK
! Input
integer nTri
real*8, dimension(nRoots*nTri) :: FMO1t
real*8, dimension(nRoots*nacpr2) :: FMO2t
! A little help from
real*8, dimension(1) :: rdum
real*8, dimension(:), allocatable :: FMO1
integer ILoc1, ILoc2, ILoc3, iOff, iS, jS, iB, jB, ji, ij, I, iK
integer ILoc4, iptmp, nConf3, L
real*8, external :: DDot_

nConf3 = nint(max(xispsm(State_SYM,1),xispsm(State_SYM,1)))
iptmp = ipGet(nconf3*nRoots)
CSFOK(:,:) = Zero
W(iptmp)%A(1:nRoots*nConf3) = Zero
call mma_allocate(FMO1,nDens)

do i=1,nRoots
  iK = I
  ILoc1 = (IK-1)*nTri
  ILoc2 = 1+(IK-1)*NACPR2
  ILoc3 = 1+(IK-1)*nConf1
  iOff = 0
  do iS=1,nSym
    jS = iS
    if (NBAS(IS)*NBAS(JS) /= 0) then
      do iB=1,nBas(iS)
        do jB=1,iB
          ioff = ioff+1
          ji = ipMat(is,js)-1+(iB-1)*nbas(iS)+jB
          FMO1(ji) = FMO1t(ioff+iLoc1)
          if (iB /= jB) then
            ij = ipMat(is,js)-1+(jB-1)*nbas(iS)+iB
            FMO1(ij) = FMO1t(ioff+iLoc1)
          end if
        end do
      end do
    end if
  end do
  !call ipin(ipCI)
  call CISigma_SA(0,State_Sym,State_Sym,FMO1,nDens,FMO2t(iLoc2),NACPR2,rdum,1,ipci,iptmp,.true.)
  CSFOK(:,IK) = CSFOK(:,IK)+real(nRoots,kind=wp)*W(iptmp)%A(iLoc3:iLoc3+nConf1-1)

  do L=1,nRoots
    ILoc4 = (L-1)*NConf1+1
    LOK(L,I) = ddot_(nConf1,CSFOK(:,IK),1,W(ipCI)%A(iLoc4),1)
  end do
end do
call mma_deallocate(FMO1)

end subroutine CalcOMat
