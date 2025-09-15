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
use MCLR_Data, only: ipCI, ipMat, nAcPr2, nConf1, nDens, XISPSM
use MCLR_procedures, only: CISigma_sa
use input_mclr, only: nBas, nRoots, nSym, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: CSFOK(nConf1,nRoots), LOK(nRoots,nRoots)
integer(kind=iwp), intent(in) :: nTri
real(kind=wp), intent(in) :: FMO1t(nTri,nRoots), FMO2t(nAcPr2,nRoots)
integer(kind=iwp) :: I, iB, ij, iOff, iptmp, iS, jB, ji, jS, L, nConf3
real(kind=wp) :: rdum(1)
real(kind=wp), allocatable :: FMO1(:)
real(kind=wp), external :: DDot_

nConf3 = nint(max(xispsm(State_SYM,1),xispsm(State_SYM,1)))
iptmp = ipGet(nconf3*nRoots)
CSFOK(:,:) = Zero
W(iptmp)%A(1:nRoots*nConf3) = Zero
call mma_allocate(FMO1,nDens)

do I=1,nRoots
  iOff = 0
  do iS=1,nSym
    jS = iS
    if (NBAS(IS)*NBAS(JS) /= 0) then
      do iB=1,nBas(iS)
        do jB=1,iB
          ioff = ioff+1
          ji = ipMat(is,js)-1+(iB-1)*nbas(iS)+jB
          FMO1(ji) = FMO1t(ioff,I)
          if (iB /= jB) then
            ij = ipMat(is,js)-1+(jB-1)*nbas(iS)+iB
            FMO1(ij) = FMO1t(ioff,I)
          end if
        end do
      end do
    end if
  end do
  !call ipin(ipCI)
  call CISigma_SA(0,State_Sym,State_Sym,FMO1,nDens,FMO2t(:,I),NACPR2,rdum,1,ipci,iptmp,.true.)
  CSFOK(:,I) = CSFOK(:,I)+real(nRoots,kind=wp)*W(iptmp)%A((I-1)*nConf1+1:I*nConf1)

  do L=1,nRoots
    LOK(L,I) = ddot_(nConf1,CSFOK(:,I),1,W(ipCI)%A((L-1)*nConf1+1),1)
  end do
end do
call mma_deallocate(FMO1)

end subroutine CalcOMat
