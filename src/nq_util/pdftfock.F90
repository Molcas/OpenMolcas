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
! Jie J. Bao, on Dec. 22, 2021, created this file.               *
! ****************************************************************
subroutine PDFTFock(FI,FA,D1,mGrid,ActMO)

use nq_pdft, only: dEdPiMO, GdEdPiMO, lft, lGGA, MOas
use nq_Info, only: mIrrep, NASHT, nIsh, nOrbt, nPot1, OffOrb
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp) :: mGrid
real(kind=wp) :: FI(nPot1), FA(nPot1), D1(NASHT**2), ActMO(mGrid*NASHT)
! Input: mGrid ActMO D1
! Output: FI FA
integer(kind=iwp) :: iGrid, iIrrep, ik, iOff1, k, nGOrb
real(kind=wp) :: dEdPiAct(mGrid*NASHT), Fact2(mGrid), SumDX(mGrid*NASHT), TempD1 !IFG
real(kind=r8), external :: DDot_

nGOrb = mGrid*nOrbt

! calculate FI. FI=2*dEdPi*pq*sum_k{kk*Fact1}
! TempD1: sum_k{kk}
! Fact2 : 2sum_k{kk}
do iGrid=1,mGrid
  TempD1 = Zero
  do iIrrep=0,mIrrep-1
    do ik=1,nIsh(iIrrep)
      k = ik+OffOrb(iIrrep)
      IOff1 = (k-1)*mGrid+iGrid
      TempD1 = TempD1+MOas(IOff1)**2
    end do
  end do
  Fact2(iGrid) = TempD1
end do

call DScal_(mGrid,Two,Fact2,1)
if (lft .and. lGGA) then
  ! In the end (drvnq_inner) there will be a process
  ! in which FI_pq=0.5(FI_pq+FI_qp)
  ! However, we do not want the factor of 0.5 for p'qrs part
  ! So we add another copy of GdEdPiMO to pick up the factor of 0.5
  call DAXpY_(nGOrb,One,GdEdPiMO,1,dEdPiMO,1)
  ! Also the pqr's part is needed with the help of the following! array
  call TransActMO2(dEdPiAct,GdEdPiMO,mGrid)
end if
call PDFTFock_Inner(FI,Fact2,dEdPiMO,MOas,mGrid)

if (lft .and. lGGA) then
  do iGrid=1,mGrid
    TempD1 = Zero
    do iIrrep=0,mIrrep-1
      do ik=1,nIsh(iIrrep)
        k = ik+OffOrb(iIrrep)
        IOff1 = (k-1)*mGrid+iGrid
        TempD1 = TempD1+MOas(IOff1)*GdEdPiMO(iOff1)
      end do
    end do
    Fact2(iGrid) = TempD1
  end do
  call DScal_(mGrid,Four,Fact2,1)
  call PDFTFock_Inner(FI,Fact2,MOas,MOas,mGrid)
end if
! calculate FA. FA=pq*sum_vx{vx*Dvx*Fact1}
! First calcualte sum_x{Dvx*x}
do iGrid=1,mGrid
  IOff1 = (iGrid-1)*NASHT+1
  call DGEMM_('T','N',NASHT,1,NASHT,One,D1,NASHT,ActMO(IOff1),NASHT,Zero,SumDX(iOff1),NASHT)
  Fact2(iGrid) = ddot_(NASHT,ActMO(iOff1),1,SumDX(iOff1),1)
end do
call PDFTFock_Inner(FA,Fact2,dEdPiMO,MOas,mGrid)
if (lft .and. lGGA) then
  do iGrid=1,mGrid
    IOff1 = (iGrid-1)*NASHT+1
    call DGEMM_('T','N',NASHT,1,NASHT,One,D1,NASHT,ActMO(IOff1),NASHT,Zero,SumDX(iOff1),NASHT)
    Fact2(iGrid) = ddot_(NASHT,dEdPiAct(iOff1),1,SumDX(iOff1),1)
  end do
  call DScal_(mGrid,Two,Fact2,1)
  call PDFTFock_Inner(FA,Fact2,MOas,MOas,mGrid)
end if

return

end subroutine PDFTFock
