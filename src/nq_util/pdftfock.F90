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
use nq_Info, only: mIrrep, NASHT, nIsh, nPot1, OffOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: FI(nPot1), FA(nPot1)
integer(kind=iwp), intent(in) :: mGrid
real(kind=wp), intent(in) :: D1(NASHT,NASHT), ActMO(NASHT,mGrid)
integer(kind=iwp) :: iGrid, iIrrep, ik, k
real(kind=wp) :: TempD1
real(kind=wp), allocatable :: dEdPiAct(:,:), Fact2(:), SumDX(:,:)
real(kind=wp), external :: DDot_

! calculate FI. FI=2*dEdPi*pq*sum_k{kk*Fact1}
! TempD1: sum_k{kk}
! Fact2 : 2sum_k{kk}
call mma_allocate(Fact2,mGrid,label='Fact2')
do iGrid=1,mGrid
  TempD1 = Zero
  do iIrrep=0,mIrrep-1
    do ik=1,nIsh(iIrrep)
      k = ik+OffOrb(iIrrep)
      TempD1 = TempD1+MOas(iGrid,k)**2
    end do
  end do
  Fact2(iGrid) = Two*TempD1
end do

if (lft .and. lGGA) then
  ! In the end (drvnq_inner) there will be a process
  ! in which FI_pq=0.5(FI_pq+FI_qp)
  ! However, we do not want the factor of 0.5 for p'qrs part
  ! So we add another copy of GdEdPiMO to pick up the factor of 0.5
  dEdPiMO(:,:) = dEdPiMO+GdEdPiMO
  ! Also the pqr's part is needed with the help of the following! array
  call mma_allocate(dEdPiAct,NASHT,mGrid,label='dEdPiAct')
  call TransActMO2(dEdPiAct,GdEdPiMO,mGrid)
end if
call PDFTFock_Inner(FI,Fact2,dEdPiMO,MOas,mGrid)

if (lft .and. lGGA) then
  do iGrid=1,mGrid
    TempD1 = Zero
    do iIrrep=0,mIrrep-1
      do ik=1,nIsh(iIrrep)
        k = ik+OffOrb(iIrrep)
        TempD1 = TempD1+MOas(iGrid,k)*GdEdPiMO(iGrid,k)
      end do
    end do
    Fact2(iGrid) = Four*TempD1
  end do
  call PDFTFock_Inner(FI,Fact2,MOas,MOas,mGrid)
end if
! calculate FA. FA=pq*sum_vx{vx*Dvx*Fact1}
! First calculate sum_x{Dvx*x}
call mma_allocate(SumDX,NASHT,mGrid,label='SumDX')
do iGrid=1,mGrid
  call DGEMM_('T','N',NASHT,1,NASHT,One,D1,NASHT,ActMO(:,iGrid),NASHT,Zero,SumDX(:,iGrid),NASHT)
  Fact2(iGrid) = ddot_(NASHT,ActMO(:,iGrid),1,SumDX(:,iGrid),1)
end do
call PDFTFock_Inner(FA,Fact2,dEdPiMO,MOas,mGrid)
if (lft .and. lGGA) then
  do iGrid=1,mGrid
    call DGEMM_('T','N',NASHT,1,NASHT,One,D1,NASHT,ActMO(:,iGrid),NASHT,Zero,SumDX(:,iGrid),NASHT)
    Fact2(iGrid) = Two*ddot_(NASHT,dEdPiAct(:,iGrid),1,SumDX(:,iGrid),1)
  end do
  call mma_deallocate(dEdPiAct)
  call PDFTFock_Inner(FA,Fact2,MOas,MOas,mGrid)
end if
call mma_deallocate(Fact2)
call mma_deallocate(SumDX)

return

end subroutine PDFTFock
