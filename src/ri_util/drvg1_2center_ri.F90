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
! Copyright (C) 1990,1991,1992,2000,2007, Roland Lindh                 *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drvg1_2Center_RI(Grad,Temp,nGrad,ij2,nij_Eff)
!***********************************************************************
!                                                                      *
!  Object: driver for 2-center two-electron integrals in the RI scheme.*
!                                                                      *
!   The integral derivative is formulated as                           *
!   -Sum(ML) X_ij^K   V_LM^(1) X_kl^L  where                           *
!                                                                      *
!  X_ij^K = Sum(L) R_ij_L  Q_L^K                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified for gradient calculation. January '92           *
!             Modified for SetUp_Ints. January '00                     *
!             Modified for 2-center RI gradients, January '07          *
!***********************************************************************

use k2_setup
use iSD_data
use pso_stuff
use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
use Basis_Info
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use RICD_Info, only: Do_RI
use Symmetry_Info, only: nIrrep
use Para_Info, only: nProcs, King
use ExTerm, only: CijK, AMP2, iMP2prpt, A

implicit real*8(A-H,O-Z)
external Rsv_Tsk
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
#include "exterm.fh"
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
integer nGrad, nij_Eff
real*8 Grad(nGrad), Temp(nGrad)
integer, allocatable :: ij2(:,:)
! Local arrays
real*8 Coor(3,4)
integer iAnga(4), iCmpa(4), iShela(4), iShlla(4), iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4)
logical EQ, Shijij, AeqB, CeqD, DoGrad, DoFock, Indexation, FreeK2, Verbose, JfGrad(3,4), ABCDeq, No_Batch, Rsv_Tsk
character format*72
integer iSD4(0:nSD,4)
save MemPrm
real*8, allocatable :: TMax2(:,:), TMax1(:), Tmp(:,:)
integer, allocatable :: Shij(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 9
iPrint = nPrint(iRout)
#ifdef _CD_TIMING_
Twoel2_CPU = 0.0d0
Twoel2_Wall = 0.0d0
Pget2_CPU = 0.0d0
Pget2_Wall = 0.0d0
#endif

iFnc(1) = 0
iFnc(2) = 0
iFnc(3) = 0
iFnc(4) = 0
PMax = Zero
Temp(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Handle only the auxiliary basis set.

if (Do_RI) then
  call Set_Basis_Mode('Auxiliary')
else
  call Set_Basis_Mode('Valence')
end if
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Precompute k2 entities.

Indexation = .false.
DoFock = .false.
DoGrad = .true.
ThrAO = Zero
call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
mSkal = nSkal
nPairs = nSkal*(nSkal+1)/2
nQuad = nPairs*(nPairs+1)/2
Pren = Zero
Prem = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
MxPrm = 0
do iAng=0,S%iAngMx
  MxPrm = max(MxPrm,S%MaxPrm(iAng))
end do
nZeta = MxPrm*MxPrm
nEta = MxPrm*MxPrm
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

if (Do_RI) then
  nTMax = nSkal
  call mma_allocate(TMax1,nTMax,Label='TMax1')
  call mma_allocate(Tmp,nSkal,nSkal,Label='Tmp')
  call Shell_MxSchwz(nSkal,Tmp)
  TMax1(1:nSkal) = Tmp(1:nSkal,nSkal)
  call mma_deallocate(Tmp)

  TMax_all = Zero
  do iS=1,nSkal-1
    TMax_all = max(TMax_all,TMax1(iS))
  end do
else
  call mma_allocate(TMax2,nSkal,nSkal,Label='TMax2')
  call Shell_MxSchwz(nSkal,TMax2)
  TMax_all = Zero
  do ij=1,nij_Eff
    iS = ij2(1,ij)
    jS = ij2(2,ij)
    TMax_all = max(TMax_all,TMax2(iS,jS))
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate some scratch arrays to be used by the pget routines.
! In particular we will have temporary arrays for A_IJ and C_ijK.

! Lower case: valence basis set
! Upper case: auxiliary basis sets

if (DoCholExch) then

  ! Find the largest number of contractions in any given shell
  ! of auxiliary functions.

  MxChVInShl = 1
  if (Do_RI) then
    do i=1,nSkal
      MxChVInShl = max(MxChVInShl,iSD(3,i))
    end do
  else
    write(6,*) 'Not Implemented for Cholesky yet!'
    call Abend()
  end if

  ! Scratch for A_IJ

  lA = MxChVInShl*MxChVInShl
  call mma_allocate(A,lA,Label='A')
  if (iMP2Prpt == 2) then
    lA_MP2 = MxChVInShl
    call mma_allocate(AMP2,lA_MP2,2,Label='AMP2')
  end if

  ! Find the largest set of ij. The basis i and j is due to the
  ! CD of the one-particle density.

  nIJRMax = 0
  do jDen=1,nKvec
    do iSym1=1,nIrrep
      do iSym2=1,nIrrep
        nIJRMax = max(nIJRMax,nIJR(iSym1,iSym2,jDen))
      end do
    end do
  end do

  ! Get scratch for C_kl^I and C_kl^J.
  ! Note that we need nDen arrays for C_kl^I and one for C_kl^J
  ! A_IJ = Sum(kl) C_kl^I x C_kl^J

  call mma_allocate(CijK,nIJRMax*MxChVInShl*(nKvec+1),Label='CijK')

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

if (Do_RI) then
  mij = (nSkal-1)
  call mma_allocate(Shij,2,mij,Label='Shij')
  nij = 0
  do iS=1,nSkal-1
    if (TMax_All*TMax1(iS) >= CutInt) then
      nij = nij+1
      Shij(1,nij) = nSkal
      Shij(2,nij) = iS
    end if
  end do
else
  mij = nij_Eff
  call mma_allocate(Shij,2,mij,Label='Shij')
  nij = 0
  do ij=1,nij_Eff
    iS = ij2(1,ij)
    jS = ij2(2,ij)
    if (TMax_All*TMax2(iS,jS) >= CutInt) then
      nij = nij+1
      Shij(1,nij) = iS
      Shij(2,nij) = jS
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute FLOP's for the transfer equation.

do iAng=0,S%iAngMx
  do jAng=0,iAng
    nHrrab = 0
    do i=0,iAng+1
      do j=0,jAng+1
        if (i+j <= iAng+jAng+1) then
          ijMax = min(iAng,jAng)+1
          nHrrab = nHrrab+ijMax*2+1
        end if
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! For a parallel implementation the iterations over shell-pairs
! are parallelized.

call Init_Tsk(id,nij*(nij+1)/2)
!                                                                      *
!***********************************************************************
!                                                                      *
! In MPP case dispatch one processor to do 1-el gradients first

if ((nProcs > 1) .and. King()) then
  if (Do_RI) call Free_iSD()
  call Drvh1(Grad,Temp,nGrad)
  !if (nPrint(1) >= 15) call PrGrad(' Gradient excluding two-electron contribution',Grad,lDisp(0),ChDisp)
  call dcopy_(nGrad,[Zero],0,Temp,1)
  if (Do_RI) then
    call Set_Basis_Mode('Auxiliary')
    call Setup_iSD()
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MemMax)
if (MemMax > 1000) MemMax = MemMax-1000
call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
ipMem1 = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
do while (Rsv_Tsk(id,jlS))
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.

  ! Now do a quadruple loop over shells

  jS_ = int((One+sqrt(Eight*dble(jlS)-Three))/Two)
  iS = Shij(1,jS_)
  jS = Shij(2,jS_)
  lS_ = int(dble(jlS)-dble(jS_)*(dble(jS_)-One)/Two)
  kS = Shij(1,lS_)
  lS = Shij(2,lS_)
  call CWTime(TCpu1,TWall1)

  if (Do_RI) then
    Aint = TMax1(jS)*TMax1(lS)
  else
    Aint = TMax2(iS,jS)*TMax2(kS,lS)
  end if
  if (AInt < CutInt) cycle
  !if ((is == 3) .and. (js == 3) .and. (ks == 1) .and. (ls == 1)) then
  !  iPrint = 15
  !  nPrint(39) = 15
  !else
  !  iPrint = nPrint(iRout)
  !  nPrint(39) = 5
  !end if
  if (iPrint >= 15) write(6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)
  call Size_SO_block_g(iSD4,nSD,nSO,No_batch)
  if (No_batch) cycle

  call Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! --------> Memory Managment <--------
  !
  ! Compute memory request for the primitives, i.e.
  ! how much memory is needed up to the transfer equation.

  call MemRys_g(iSD4,nSD,nRys,MemPrm)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ABCDeq = EQ(Coor(1,1),Coor(1,2)) .and. EQ(Coor(1,1),Coor(1,3)) .and. EQ(Coor(1,1),Coor(1,4))
  ijklA = iSD4(1,1)+iSD4(1,2)+iSD4(1,3)+iSD4(1,4)
  if ((nIrrep == 1) .and. ABCDeq .and. (mod(ijklA,2) == 1)) cycle
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Decide on the partioning of the shells based on the
  ! available memory and the requested memory.
  !
  ! Now check if all blocks can be computed and stored at once.

  call SOAO_g(iSD4,nSD,nSO,MemPrm,MemMax,iBsInc,jBsInc,kBsInc,lBsInc,iPrInc,jPrInc,kPrInc,lPrInc,ipMem1,ipMem2,Mem1,Mem2,iFnc, &
              MemPSO)
  iBasi = iSD4(3,1)
  jBasj = iSD4(3,2)
  kBask = iSD4(3,3)
  lBasl = iSD4(3,4)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Int_Parm_g(iSD4,nSD,iAnga,iCmpa,iShlla,iShela,iPrimi,jPrimj,kPrimk,lPriml,k2ij,nDCRR,k2kl,nDCRS,mdci,mdcj,mdck,mdcl,AeqB, &
                  CeqD,nZeta,nEta,ipZeta,ipZI,ipP,ipEta,ipEI,ipQ,ipiZet,ipiEta,ipxA,ipxB,ipxG,ipxD,l2DI,nab,nHmab,ncd,nHmcd,nIrrep)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scramble arrays (follow angular index)

  do iCar=1,3
    do iSh=1,4
      JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
      if (((iSh == 1) .or. (iSh == 3)) .and. Do_RI) then
        JfGrad(iCar,iSh) = .false.
      else if (iand(iSD4(15,iSh),2**(iCar-1)) == 2**(iCar-1)) then
        JfGrad(iCar,iSh) = .true.
      else
        JfGrad(iCar,iSh) = .false.
      end if
    end do
  end do

  do iBasAO=1,iBasi,iBsInc
    iBasn = min(iBsInc,iBasi-iBasAO+1)
    iAOst(1) = iBasAO-1
    do jBasAO=1,jBasj,jBsInc
      jBasn = min(jBsInc,jBasj-jBasAO+1)
      iAOst(2) = jBasAO-1

      do kBasAO=1,kBask,kBsInc
        kBasn = min(kBsInc,kBask-kBasAO+1)
        iAOst(3) = kBasAO-1
        do lBasAO=1,lBasl,lBsInc
          lBasn = min(lBsInc,lBasl-lBasAO+1)
          iAOst(4) = lBasAO-1

          ! Get the 2nd order density matrix in SO basis.

          nijkl = iBasn*jBasn*kBasn*lBasn
#         ifdef _CD_TIMING_
          call CWTIME(Pget0CPU1,Pget0WALL1)
#         endif
          call PGet0(iCmpa,iBasn,jBasn,kBasn,lBasn,Shijij,iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,iFnc(1)*iBasn,iFnc(2)*jBasn, &
                     iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
#         ifdef _CD_TIMING_
          call CWTIME(Pget0CPU2,Pget0WALL2)
          Pget2_CPU = Pget2_CPU+Pget0CPU2-Pget0CPU1
          Pget2_Wall = Pget2_Wall+Pget0WALL2-Pget0WALL1
#         endif
          if (AInt*PMax < CutInt) cycle

          ! Compute gradients of shell quadruplet

#         ifdef _CD_TIMING_
          call CWTIME(TwoelCPU1,TwoelWall1)
#         endif
          call TwoEl_g(Coor,iAnga,iCmpa,iShela,iShlla,iAOV,mdci,mdcj,mdck,mdcl,nRys,Data_k2(k2ij),nab,nHmab,nDCRR,Data_k2(k2kl), &
                       ncd,nHmcd,nDCRS,Pren,Prem,iPrimi,iPrInc,jPrimj,jPrInc,kPrimk,kPrInc,lPriml,lPrInc, &
                       Shells(iSD4(0,1))%pCff(1,iBasAO),iBasn,Shells(iSD4(0,2))%pCff(1,jBasAO),jBasn, &
                       Shells(iSD4(0,3))%pCff(1,kBasAO),kBasn,Shells(iSD4(0,4))%pCff(1,lBasAO),lBasn,Mem_DBLE(ipZeta), &
                       Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,Mem_DBLE(ipEta),Mem_DBLE(ipEI),Mem_DBLE(ipQ),nEta,Mem_DBLE(ipxA), &
                       Mem_DBLE(ipxB),Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,JfGrad,JndGrd,Sew_Scr(ipMem1),nSO,Sew_Scr(ipMem2), &
                       Mem2,Aux,nAux,Shijij)
#         ifdef _CD_TIMING_
          call CWTIME(TwoelCPU2,TwoelWall2)
          Twoel2_CPU = Twoel2_CPU+TwoelCPU2-TwoelCPU1
          Twoel2_Wall = Twoel2_Wall+TwoelWall2-TwoelWall1
#         endif
          if (iPrint >= 15) call PrGrad(' In Drvg1_2Center_RI: Grad',Temp,nGrad,ChDisp)

        end do
      end do

    end do
  end do

end do
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Sew_Scr)
call Free_Tsk(id)
call mma_deallocate(Shij)
if (allocated(TMax1)) call mma_deallocate(TMax1)
if (allocated(TMax2)) call mma_deallocate(TMax2)
!                                                                      *
!***********************************************************************
!                                                                      *
Verbose = .false.
FreeK2 = .true.
call Term_Ints(Verbose,FreeK2)
!                                                                      *
!***********************************************************************
!                                                                      *
call Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)

iPren = 3+max(1,int(log10(Pren+0.001D+00)))
iPrem = 3+max(1,int(log10(Prem+0.001D+00)))
write(format,'(A,I2,A,I2,A)') '(A,F',iPren,'.0,A,F',iPrem,'.0,A)'
if (iPrint >= 6) then
  write(6,format) ' A total of',Pren,' entities were prescreened and',Prem,' were kept.'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (DoCholExch) then
  call mma_deallocate(CijK)
  call mma_deallocate(A)
end if
if (allocated(AMP2)) call mma_deallocate(AMP2)

call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Drvg1_2Center_RI
