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

use setup, only: mSkal, MxPrm, nAux
use Index_Functions, only: nTri_Elem
use iSD_data, only: iSD, nSD
use pso_stuff, only: A_PT2, nBasA, nBasASQ, nBasT
use k2_arrays, only: Aux, Destroy_BraKet, Sew_Scr
use k2_structure, only: k2data
use Disp, only: ChDisp, l2DI
use Basis_Info, only: nBas, nBas_Aux, Shells
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use RICD_Info, only: Do_RI
use Symmetry_Info, only: nIrrep
use Para_Info, only: King, nProcs
use RI_glob, only: A, AMP2, CijK, DoCholExch, iMP2prpt, MxChVInShl, nIJR, nKvec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Eight, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad, nij_Eff, ij2(2,nij_Eff)
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
#include "print.fh"
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
integer(kind=iwp) :: i, iAng, iAnga(4), iAOst(4), iAOV(4), iBasAO, iBasi, iBasn, iBsInc, iCar, iCmpa(4), id, iFnc(4), ij, ijkla, &
                     ijMax, ik2, ipMem1, ipMem2, iPrem, iPren, iPrimi, iPrInc, iPrint, iRout, iS, iSD4(0:nSD,4), iSh, iShela(4), &
                     iShlla(4), istabs(4), iSym1, iSym2, j, jAng, jBasAO, jBasj, jBasn, jBsInc, jDen, jk2, jlS, JndGrd(3,4), &
                     jPrimj, jPrInc, jS, jS_, k2ij, k2kl, kBasAO, kBask, kBasn, kBsInc, kBtch, kPrimk, kPrInc, kS, lA, lA_MP2, &
                     lBasAO, lBasl, lBasn, lBsInc, lPriml, lPrInc, lS, lS_, LUAPT2, mBtch, mdci, mdcj, mdck, mdcl, Mem1, Mem2, &
                     MemMax, MemPSO, mij, nab, nBtch, ncd, nDCRR, nDCRS, nEta, nHmab, nHMcd, nHrrab, nij, nijkl, nIJRMax, nPairs, &
                     nQuad, nRys, nSkal, nSO, nTMax, nZeta
real(kind=wp) :: A_int, Coor(3,4), PMax, Prem, Pren, TCpu1, ThrAO, TMax_all, TWall1
#ifdef _CD_TIMING_
real(kind=wp) :: Pget0CPU1, Pget0CPU2, Pget0WALL1, Pget0WALL2, TwoelCPU1, TwoelCPU2, TwoelWall1, TwoelWall2
#endif
logical(kind=iwp) :: ABCDeq, AeqB, CeqD, DoFock, DoGrad, EQ, Indexation, JfGrad(3,4), No_Batch, Shijij
character(len=72) :: frmt
character(len=8) :: Method_chk
integer(kind=iwp), save :: MemPrm
integer(kind=iwp), allocatable :: Shij(:,:)
real(kind=wp), allocatable :: TMax1(:), TMax2(:,:), Tmp(:,:)
integer(kind=iwp), external :: IsFreeUnit
logical(kind=iwp), external :: Rsv_Tsk

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 9
iPrint = nPrint(iRout)
#ifdef _CD_TIMING_
Twoel2_CPU = Zero
Twoel2_Wall = Zero
Pget2_CPU = Zero
Pget2_Wall = Zero
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
nPairs = nTri_Elem(nSkal)
nQuad = nTri_Elem(nPairs)
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
    write(u6,*) 'Not Implemented for Cholesky yet!'
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
  mij = nSkal-1
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
! CASPT2

call Get_cArray('Relax Method',Method_chk,8)
if (Method_chk == 'CASPT2  ') then
  ! Just read A_{JK} type matrix constructed in CASPT2
  nBasT = 0
  nBasA = 0
  nBasASQ = 0
  do iSym1=0,nIrrep-1
    nBasT = nBasT+nBas(iSym1)
    nBasA = nBasA+nBas_Aux(iSym1)-1
    nBasASQ = nBasASQ+(nBas_Aux(iSym1)-1)**2
  end do
  call mma_allocate(A_PT2,nBasA,nBasA,Label='A_PT2')
  ! Now, read
  !call PrgmTranslate('CMOPT2',RealName,lRealName)
  !LuCMOPT2 = 61
  !call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.false.,1,'OLD',is_error)
  !read(LuCMOPT2) A_PT2(1:nBasASQ,1)
  !close(LuCMOPT2)

  ! Read A_PT2 from LUAPT2
  LuAPT2 = isFreeUnit(68)
  call daname_mf_wa(LUAPT2,'A_PT2')
  id = 0
  call ddafile(LUAPT2,2,A_PT2,nBasASq,id)
  call daclos(LUAPT2)
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

call Init_Tsk(id,nTri_Elem(nij))
!                                                                      *
!***********************************************************************
!                                                                      *
! In MPP case dispatch one processor to do 1-el gradients first

if ((nProcs > 1) .and. King()) then
  if (Do_RI) call Free_iSD()
  call Drvh1(Grad,Temp,nGrad)
  !if (nPrint(1) >= 15) call PrGrad(' Gradient excluding two-electron contribution',Grad,lDisp(0),ChDisp)
  Temp(:) = Zero
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

  jS_ = int((One+sqrt(Eight*real(jlS,kind=wp)-Three))*Half)
  iS = Shij(1,jS_)
  jS = Shij(2,jS_)
  lS_ = int(real(jlS,kind=wp)-real(jS_,kind=wp)*(real(jS_,kind=wp)-One)*Half)
  kS = Shij(1,lS_)
  lS = Shij(2,lS_)
  call CWTime(TCpu1,TWall1)

  if (Do_RI) then
    A_int = TMax1(jS)*TMax1(lS)
  else
    A_int = TMax2(iS,jS)*TMax2(kS,lS)
  end if
  if (A_Int < CutInt) cycle
  !if ((is == 3) .and. (js == 3) .and. (ks == 1) .and. (ls == 1)) then
  !  iPrint = 15
  !  nPrint(39) = 15
  !else
  !  iPrint = nPrint(iRout)
  !  nPrint(39) = 5
  !end if
  if (iPrint >= 15) write(u6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
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
  call Int_Parm_g(iSD4,nSD,iAnga,iCmpa,iShlla,iShela,iPrimi,jPrimj,kPrimk,lPriml, &
                  k2ij,ik2,nDCRR,k2kl,jk2,nDCRS,mdci,mdcj,mdck,mdcl, &
                  AeqB,CeqD,nZeta,nEta,l2DI,nab,nHmab,ncd,nHmcd,nIrrep)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scramble arrays (follow angular index)

  do iCar=1,3
    do iSh=1,4
      JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
      if (((iSh == 1) .or. (iSh == 3)) .and. Do_RI) then
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
          if (A_Int*PMax < CutInt) cycle

          ! Compute gradients of shell quadruplet

#         ifdef _CD_TIMING_
          call CWTIME(TwoelCPU1,TwoelWall1)
#         endif
          call TwoEl_g(Coor,iAnga,iCmpa,iShela,iShlla,iAOV,mdci,mdcj,mdck,mdcl,nRys, &
                       k2Data(:,ik2),k2Data(:,jk2), &
                       nDCRR,nDCRS,Pren,Prem,iPrimi,iPrInc,jPrimj,jPrInc,kPrimk,kPrInc,lPriml,lPrInc, &
                       Shells(iSD4(0,1))%pCff(1,iBasAO),iBasn,Shells(iSD4(0,2))%pCff(1,jBasAO),jBasn, &
                       Shells(iSD4(0,3))%pCff(1,kBasAO),kBasn,Shells(iSD4(0,4))%pCff(1,lBasAO),lBasn, &
                       nZeta,nEta,Temp,nGrad,JfGrad,JndGrd,Sew_Scr(ipMem1),nSO, &
                       Sew_Scr(ipMem2),Mem2,Aux,nAux,Shijij)
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

  call Destroy_Braket()

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
if (Method_chk == 'CASPT2') call mma_deallocate(A_PT2)
call mma_deallocate(Shij)
if (allocated(TMax1)) call mma_deallocate(TMax1)
if (allocated(TMax2)) call mma_deallocate(TMax2)
!                                                                      *
!***********************************************************************
!                                                                      *
call Term_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
call Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)

iPren = 3+max(1,int(log10(Pren+0.001_wp)))
iPrem = 3+max(1,int(log10(Prem+0.001_wp)))
write(frmt,'(A,I2,A,I2,A)') '(A,F',iPren,'.0,A,F',iPrem,'.0,A)'
if (iPrint >= 6) then
  write(u6,frmt) ' A total of',Pren,' entities were prescreened and',Prem,' were kept.'
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
