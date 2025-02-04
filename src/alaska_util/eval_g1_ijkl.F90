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
! Copyright (C) 1990-1992,2000,2025, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

subroutine Eval_g1_ijkl(iS,jS,kS,lS,Temp,nGrad,A_Int)

use Index_Functions, only: nTri_Elem
use k2_arrays, only: Create_BraKet, Destroy_BraKet, Sew_Scr
use iSD_data, only: iSD, nSD
use Symmetry_Info, only: nIrrep
use Gateway_Info, only: CutInt
use setup, only: mSkal
use RICD_Info, only: RI_2C, RI_3C
use stdalloc, only: mma_allocate, mma_maxDBLE
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iS, jS, kS, lS, nGrad
real(kind=wp), intent(inout) :: Temp(nGrad)
real(kind=wp), intent(in) :: A_Int
integer(kind=iwp) :: iBasAO, iBasi, iBasn, iBsInc, iCar, iFnc(4), ijklA, ipMem1, ipMem2, iSD4(0:nSD,4), iSh, jBasAO, jBasj, jBasn, &
                     jBsInc, JndGrd(3,4), kBasAO, kBask, kBasn, kBsInc, lBasAO, lBasl, lBasn, lBsInc, Mem1, Mem2, MemMax, MemPrm, &
                     MemPSO, nijkl, nPairs, nQuad, nRys, nSO
real(kind=wp) :: Coor(3,4), PMax
real(kind=wp), pointer :: PSO(:), Wrk2(:)
logical(kind=iwp) :: ABCDeq, JfGrad(3,4), No_batch
integer(kind=iwp), external :: MemSO2_P
logical(kind=iwp), external :: EQ

iFnc(:) = 0
PMax = Zero
nPairs = nTri_Elem(mSkal)
nQuad = nTri_Elem(nPairs)

if (.not. allocated(Sew_Scr)) then
  call mma_MaxDBLE(MemMax)
  if (MemMax > 8000) MemMax = MemMax-8000
  call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
else
  MemMax = size(Sew_Scr)
end if
ipMem1 = 1

#ifdef _DEBUGPRINT_
write(u6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
#endif
!                                                                *
!*****************************************************************
!                                                                *
call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)

nSO = MemSO2_P(nSD,iSD4)
No_batch = (nSO == 0)
if (No_batch) return

call Coor_Setup(iSD4,nSD,Coor)
!                                                                      *
!***********************************************************************
!                                                                      *
! -------> Memory Managment <--------
!
! Compute memory request for the primitives, i.e.
! how much memory is needed up to the transfer
! equation.

call MemRys_g(iSD4,nSD,nRys,MemPrm)
!                                                                      *
!***********************************************************************
!                                                                      *
ABCDeq = EQ(Coor(1,1),Coor(1,2)) .and. EQ(Coor(1,1),Coor(1,3)) .and. EQ(Coor(1,1),Coor(1,4))
ijklA = iSD4(1,1)+iSD4(1,2)+iSD4(1,3)+iSD4(1,4)
if ((nIrrep == 1) .and. ABCDeq .and. (mod(ijklA,2) == 1)) return

!                                                                      *
!***********************************************************************
!                                                                      *
! partition memory for K2(ij)/K2(kl) temp spaces zeta,eta,kappa,P,Q

call Create_BraKet(iSD4(5,1)*iSD4(5,2),iSD4(5,3)*iSD4(5,4))

!                                                                      *
!***********************************************************************
!                                                                      *
! Decide on the partitioning of the shells based on the
! available memory and the requested memory.
!
! Now check if all blocks can be computed and stored at once.

call PSOAO1(nSO,MemPrm,MemMax,iFnc,ipMem1,ipMem2,Mem1,Mem2,MemPSO,nSD,iSD4)

PSO(1:Mem1) => Sew_Scr(ipMem1:ipMem1+Mem1-1)
Wrk2(1:Mem2) => Sew_Scr(ipMem2:ipMem2+Mem2-1)

iBasi = iSD4(3,1)
jBasj = iSD4(3,2)
kBask = iSD4(3,3)
lBasl = iSD4(3,4)

iBsInc = iSD4(4,1)
jBsInc = iSD4(4,2)
kBsInc = iSD4(4,3)
lBsInc = iSD4(4,4)

!                                                                      *
!***********************************************************************
!                                                                      *
! Scramble arrays (follow angular index)

do iCar=1,3
  do iSh=1,4
    JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
    if (RI_3C .and. (iSh == 1)) then
      JfGrad(iCar,iSh) = .false.
      JndGrd(iCar,iSh) = 0
    else if (RI_2C .and. ((iSh == 1) .or. (iSh == 3))) then
      JfGrad(iCar,iSh) = .false.
    else if (btest(iSD4(15,iSh),iCar-1)) then
      JfGrad(iCar,iSh) = .true.
    else
      JfGrad(iCar,iSh) = .false.
    end if
  end do
end do

do iBasAO=1,iBasi,iBsInc
  iBasn = min(iBsInc,iBasi-iBasAO+1)
  iSD4(8,1) = iBasAO-1
  iSD4(19,1) = iBasn

  do jBasAO=1,jBasj,jBsInc
    jBasn = min(jBsInc,jBasj-jBasAO+1)
    iSD4(8,2) = jBasAO-1
    iSD4(19,2) = jBasn

    do kBasAO=1,kBask,kBsInc
      kBasn = min(kBsInc,kBask-kBasAO+1)
      iSD4(8,3) = kBasAO-1
      iSD4(19,3) = kBasn

      do lBasAO=1,lBasl,lBsInc
        lBasn = min(lBsInc,lBasl-lBasAO+1)
        iSD4(8,4) = lBasAO-1
        iSD4(19,4) = lBasn

        ! Get the 2nd order density matrix in SO basis.

        nijkl = iBasn*jBasn*kBasn*lBasn

        ! Fetch the T_i,j,kappa, lambda corresponding to
        ! kappa = k, lambda = l

        call PGet0(nijkl,PSO,nSO,iFnc,MemPSO,Wrk2,Mem2,nQuad,PMax,iSD4)
        if (A_Int*PMax < CutInt) return

        ! Compute gradients of shell quadruplet

        call TwoEl_g(Coor,nRys,Temp,nGrad,JfGrad,JndGrd,PSO,nijkl,nSO,Sew_Scr(ipMem2),Mem2,iSD4)

#       ifdef _DEBUGPRINT_
        call PrGrad(' In Drvg1: Grad',Temp,nGrad)
#       endif

      end do
    end do

  end do
end do

call Destroy_BraKet()

end subroutine Eval_g1_ijkl
