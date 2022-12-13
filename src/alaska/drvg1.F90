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
! Copyright (C) 1990-1992,2000, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drvg1(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals. The four outermost loops *
!          will control the type of the two-electron integral, eg.     *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis        *
!          functions of the requested type.                            *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified for gradient calculation. January '92           *
!             Modified for SetUp_Ints. January '00                     *
!***********************************************************************

use k2_setup, only: Data_k2
use iSD_data, only: iSD
use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
use Aces_Stuff, only: G_toc, nSSDM, SSDM
use Basis_Info, only: nBas, Shells
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Para_Info, only: nProcs, King
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Eight
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
#include "Molcas.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
integer(kind=iwp) :: i, iAng, iAnga(4), iAOst(4), iAOV(4), iBasAO, iBasi, iBasn, iBsInc, iCar, iCmpa(4), iCnt, iFnc(4), ijklA, &
                     ijMax, ijS, iOpt, iost, ipEI, ipiEta, ipMem1, ipMem2, ipP, ipQ, iPrem, iPren, iPrimi, iPrInc, iPrint, ipEta, &
                     ipxA, ipxB, ipxD, ipxG, ipZI, iRout, iS, iSD4(0:nSD,4), iSh, iShela(4), iShlla(4), iSSDM, istabs(4), j, jAng, &
                     jBAsAO, jBasj, jBasn, jBsInc, jPrimj, jPrInc, jS, JndGrd(3,4), k2ij, k2kl, kBasAO, kBask, kBasn, kBsInc, &
                     kBtch, kls, kPrimk, kPrInc, kS, lBasAO, lBasl, lBasn, lBsInc, lPriml, lPrInc, lRealName, lS, luGamma, &
                     luCMOPT2, MaxShlAO, mBtch, mdci, mdcj, mdck, mdcl, Mem1, Mem2, MemMax, MemPSO, nab, nBasI, nBasT, nBtch, ncd, &
                     nDCRR, nDCRS, nEta, nFro(8), nHmab, nHmcd, nHrrab, nij, nijkl, nOcc(8), nPairs, nQuad, nRys, nSkal, nSO, nZeta
real(kind=wp) :: A_int, Cnt, Coor(3,4), P_Eff, PMax, Prem, Pren, TCpu1, TCpu2, ThrAO, TMax_all, TskHi, TskLw, TWall1, TWall2
logical(kind=iwp) :: ABCDeq, AeqB, CeqD, DoFock, DoGrad, EQ, FreeK2, Indexation, is_error, JfGrad(3,4), lDummy, Loadvec, No_Batch, &
                     Shijij, Skip, Triangular, Verbose
character(len=4096) :: RealName
character(len=72) :: formt
character(len=8) :: Method_chk
integer(kind=iwp), allocatable :: Ind_ij(:,:), iOffAO(:)
real(kind=wp), allocatable :: CMOPT2(:), TMax(:,:), WRK1(:), WRK2(:)
integer(kind=iwp), save :: MemPrm
logical(kind=iwp), external :: Rsv_GTList
!*********** columbus interface ****************************************
integer(kind=iwp) :: Columbus
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
real(kind=wp) :: Pget0CPU1, Pget0CPU2, Pget0WALL1, Pget0WALL2, Pget_CPU, Pget_Wall, Total_Dens_CPU, Total_Dens_Wall, &
                 Total_Der_CPU, Total_Der_CPU2, Total_Der_Wall, Total_Der_Wall2, Twoel_CPU, Twoel_Wall, TwoelCPU1, TwoelCPU2, &
                 TwoelWall1, TwoelWall2
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 9
iPrint = nPrint(iRout)
iPrint = 0

iFnc(1) = 0
iFnc(2) = 0
iFnc(3) = 0
iFnc(4) = 0
PMax = Zero
#ifdef _CD_TIMING_
Twoel_CPU = Zero
Twoel_Wall = Zero
Pget_CPU = Zero
Pget_Wall = Zero
#endif
Temp(:) = Zero

call StatusLine(' Alaska:',' Computing 2-electron gradients')
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
!-- Precompute k2 entities.

call Get_iScalar('Columbus',Columbus)
Indexation = .false.
! MP2 gradients:
call Get_cArray('Relax Method',Method_chk,8)
if (Method_chk == 'MBPT2   ') Indexation = .true.
!*********** columbus interface ****************************************
! in order to access the half-sorted density matrix file
! some index arrays must be generated
if (Columbus == 1) Indexation = .true.
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
!-- Prepare handling of two-particle density.

call PrepP

if (Method_chk == 'CASPT2  ') then
  nBasT = 0
  do i=0,nIrrep-1
    nBasT = nBasT+nBas(i)
  end do
  nSSDM = 0

  !! The two MO indices in the half-transformed amplitude are
  !! not CASSCF but quasi-canonical orbitals.
  call mma_allocate(CMOPT2,nBasT*nBasT,Label='CMOPT2')
  call PrgmTranslate('CMOPT2',RealName,lRealName)
  LuCMOPT2 = 61
  call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.false.,1,'OLD',is_error)
  do i=1,nBasT*nBasT
    read(LuCMOPT2) CMOPT2(i)
  end do
  read(LuCMOPT2) nOcc(1)
  read(LuCMOPT2) nOcc(2)
  read(LuCMOPT2) nOcc(3)
  read(LuCMOPT2) nOcc(4)
  read(LuCMOPT2) nOcc(5)
  read(LuCMOPT2) nOcc(6)
  read(LuCMOPT2) nOcc(7)
  read(LuCMOPT2) nOcc(8)
  read(LuCMOPT2) nFro(1)
  read(LuCMOPT2) nFro(2)
  read(LuCMOPT2) nFro(3)
  read(LuCMOPT2) nFro(4)
  read(LuCMOPT2) nFro(5)
  read(LuCMOPT2) nFro(6)
  read(LuCMOPT2) nFro(7)
  read(LuCMOPT2) nFro(8)
  read(LuCMOPT2) nSSDM

  if (nSSDM /= 0) then
    call mma_allocate(SSDM,nBas(0)*(nBas(0)+1)/2,2,nSSDM,Label='SSDM')
    do iSSDM=1,nSSDM
      do i=1,nBas(0)*(nBas(0)+1)/2
        read(LuCMOPT2) SSDM(i,1,iSSDM),SSDM(i,2,iSSDM)
      end do
    end do
  end if
  close(LuCMOPT2)
  write(u6,*) 'Number of Non-Frozen Occupied Orbitals = ',nOcc(1)
  write(u6,*) 'Number of     Frozen          Orbitals = ',nFro(1)

  call mma_allocate(iOffAO,nSkal+1,Label='iOffAO')
  MaxShlAO = 0
  iOffAO(1) = 0
  do iSh=1,nSkal
    nBasI = iSD(2,iSh)*iSD(3,iSh)
    if (nBasI > MaxShlAO) MaxShlAO = nBasI
    iOffAO(iSh+1) = iOffAO(iSh)+nBasI
  end do
  call mma_allocate(G_toc,MaxShlAO**4,Label='GtocCASPT2')

  call PrgmTranslate('GAMMA',RealName,lRealName)
  LuGamma = 60
  call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nOcc(1)*nOcc(1)*8,'OLD',is_error)

  call mma_allocate(WRK1,nOcc(1)*nOcc(1),Label='WRK1')
  call mma_allocate(WRK2,MaxShlAO*nOcc(1),Label='WRK2')
end if
!                                                                      *
!***********************************************************************
!                                                                      *
MxPrm = 0
do iAng=0,S%iAngMx
  MxPrm = max(MxPrm,S%MaxPrm(iAng))
end do
nZeta = MxPrm*MxPrm
nEta = MxPrm*MxPrm
!
!***********************************************************************
!                                                                      *
!-- Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
call Shell_MxSchwz(nSkal,TMax)
TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(Ind_ij,2,nskal*(nSkal+1)/2,Label='Ind_ij')
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if ((TMax_All*TMax(iS,jS) >= CutInt) .or. (Method_chk == 'CASPT2  ')) then
      nij = nij+1
      Ind_ij(1,nij) = iS
      Ind_ij(2,nij) = jS
    end if
  end do
end do
P_Eff = real(nij,kind=wp)
!                                                                      *
!***********************************************************************
!                                                                      *
!-- Compute FLOPs for the transfer equation.

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
Triangular = .true.
call Init_TList(Triangular,P_Eff)
call Init_PPList
call Init_GTList
iOpt = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! In MPP case dispatch one processor to do 1-el gradients first

if ((nProcs > 1) .and. King()) then
  call Drvh1(Grad,Temp,nGrad)
  !if (nPrint(1) >= 15) call PrGrad(' Gradient excluding two-electron contribution',Grad,lDisp(0),ChDisp)
  Temp(:) = Zero
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MemMax)
if (MemMax > 1000) MemMax = MemMax-1000
call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
ipMem1 = 1
ijS = 0
klS = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
do
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.
  if (.not. Rsv_GTList(TskLw,TskHi,iOpt,lDummy)) exit

  ! Now do a quadruple loop over shells

  call Get_cArray('Relax Method',Method_chk,8)
  if (Method_chk /= 'CASPT2  ') then
    ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
    iS = Ind_ij(1,ijS)
    jS = Ind_ij(2,ijS)
    klS = int(TskLw-real(ijS,kind=wp)*(real(ijS,kind=wp)-One)/Two)
    kS = Ind_ij(1,klS)
    lS = Ind_ij(2,klS)
  else
    iS = 1
    jS = 1
    kS = 1
    lS = 1
    !! proceed the index
    do iCnt=1,int(TskLw)-1
      call CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)
    end do
    Cnt = real(iCnt,kind=wp)
    !! If LoadVec is true, a new vector of the half-transformed
    !! T-amplitude is read. In the first loop, it is always true.
    !! In other loops, a new vector is read only when I- and K-th
    !! are different from the previous loop.
    !! The half back-transformation, T_{ij}^{ab} ->
    !! T_{ij}^{rho sigma}, is done somewhere in CASPT2.
    !! rho and sigma correspond to either I- or K-th shells.
    !! Occupied orbital indices (correspond to J- or L-th shells)
    !! are back-transformed on-the-fly.
    LoadVec = .true.
  end if
  Cnt = TskLw
  call CWTime(TCpu1,TWall1)

  do
    A_int = TMax(iS,jS)*TMax(kS,lS)
    Skip = .false.
    if (A_Int < CutInt) Skip = .true.
    if (.not. Skip) then
      if (iPrint >= 15) write(u6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
      !                                                                *
      !*****************************************************************
      !                                                                *
      call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)
      call Size_SO_block_g(iSD4,nSD,nSO,No_batch)
      if (No_batch) Skip = .true.
    end if

    if (.not. Skip) then
      call Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! -------> Memory Managment <--------
      !
      ! Compute memory request for the primitives, i.e.
      ! how much memory is needed up to the transfer
      ! equation.

      call MemRys_g(iSD4,nSD,nRys,MemPrm)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ABCDeq = EQ(Coor(1,1),Coor(1,2)) .and. EQ(Coor(1,1),Coor(1,3)) .and. EQ(Coor(1,1),Coor(1,4))
      ijklA = iSD4(1,1)+iSD4(1,2)+iSD4(1,3)+iSD4(1,4)
      if ((nIrrep == 1) .and. ABCDeq .and. (mod(ijklA,2) == 1)) Skip = .true.
    end if
    if (.not. Skip) then
      !                                                                *
      !*****************************************************************
      !                                                                *
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
      !                                                                *
      !*****************************************************************
      !                                                                *
      call Int_Parm_g(iSD4,nSD,iAnga,iCmpa,iShlla,iShela,iPrimi,jPrimj,kPrimk,lPriml,k2ij,nDCRR,k2kl,nDCRS,mdci,mdcj,mdck,mdcl, &
                      AeqB,CeqD,nZeta,nEta,ipZeta,ipZI,ipP,ipEta,ipEI,ipQ,ipiZet,ipiEta,ipxA,ipxB,ipxG,ipxD,l2DI,nab,nHmab,ncd, &
                      nHmcd,nIrrep)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Scramble arrays (follow angular index)

      do iCar=1,3
        do iSh=1,4
          JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
          if (btest(iSD4(15,iSh),iCar-1)) then
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

              ! Fetch the T_i,j,kappa, lambda corresponding to
              ! kappa = k, lambda = l

#             ifdef _CD_TIMING_
              call CWTIME(Pget0CPU1,Pget0WALL1)
#             endif
              if (Method_chk == 'CASPT2  ') call CASPT2_BTAMP(iS,jS,kS,lS,iFnc(1)*iBasn,iFnc(2)*jBasn,iFnc(3)*kBasn,iFnc(4)*lBasn, &
                                                              iOffAO,nBasT,nOcc(1),CMOPT2(1+nbast*nfro(1)),WRK1,WRK2,G_Toc)
              call PGet0(iCmpa,iBasn,jBasn,kBasn,lBasn,Shijij,iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,iFnc(1)*iBasn,iFnc(2)*jBasn, &
                         iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
              if (A_Int*PMax < CutInt) cycle
#             ifdef _CD_TIMING_
              call CWTIME(Pget0CPU2,Pget0WALL2)
              Pget_CPU = Pget_CPU+Pget0CPU2-Pget0CPU1
              Pget_Wall = Pget_Wall+Pget0WALL2-Pget0WALL1
#             endif
              if (A_Int*PMax < CutInt) cycle

              ! Compute gradients of shell quadruplet

#             ifdef _CD_TIMING_
              call CWTIME(TwoelCPU1,TwoelWall1) ! timing_cdscf
#             endif
              call TwoEl_g(Coor,iAnga,iCmpa,iShela,iShlla,iAOV,mdci,mdcj,mdck,mdcl,nRys,Data_k2(k2ij),nab,nHmab,nDCRR, &
                           Data_k2(k2kl),ncd,nHmcd,nDCRS,Pren,Prem,iPrimi,iPrInc,jPrimj,jPrInc,kPrimk,kPrInc,lPriml,lPrInc, &
                           Shells(iSD4(0,1))%pCff(1,iBasAO),iBasn,Shells(iSD4(0,2))%pCff(1,jBasAO),jBasn, &
                           Shells(iSD4(0,3))%pCff(1,kBasAO),kBasn,Shells(iSD4(0,4))%pCff(1,lBasAO),lBasn,Mem_DBLE(ipZeta), &
                           Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,Mem_DBLE(ipEta),Mem_DBLE(ipEI),Mem_DBLE(ipQ),nEta,Mem_DBLE(ipxA), &
                           Mem_DBLE(ipxB),Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,JfGrad,JndGrd,Sew_Scr(ipMem1),nSO, &
                           Sew_Scr(ipMem2),Mem2,Aux,nAux,Shijij)
#             ifdef _CD_TIMING_
              call CWTIME(TwoelCPU2,TwoelWall2)
              Twoel_CPU = Twoel_CPU+TwoelCPU2-TwoelCPU1
              Twoel_Wall = Twoel_Wall+TwoelWall2-TwoelWall1
#             endif
              if (iPrint >= 15) call PrGrad(' In Drvg1: Grad',Temp,nGrad,ChDisp)

            end do
          end do

        end do
      end do

    end if

    Cnt = Cnt+One
    if (Cnt-TskHi > 1.0e-10_wp) then
      exit
    else if (Method_chk == 'CASPT2  ') then
      call CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)
    else
      klS = klS+1
      if (klS > ijS) then
        ijS = ijS+1
        klS = 1
      end if
      iS = Ind_ij(1,ijS)
      jS = Ind_ij(2,ijS)
      kS = Ind_ij(1,klS)
      lS = Ind_ij(2,klS)
    end if
  end do

  ! Task endpoint
  call CWTime(TCpu2,TWall2)
end do
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _MOLCAS_MPP_
call GADGOP(Temp,nGrad,'+')
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Sew_Scr)
call Free_GTList
call Free_PPList
call Free_TList
call mma_deallocate(Ind_ij)
call mma_deallocate(TMax)
if (Method_chk == 'CASPT2  ') then
  close(LuGamma)
  call mma_deallocate(iOffAO)
  call mma_deallocate(CMOPT2)
  if (nSSDM /= 0) call mma_deallocate(SSDM)
  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call CloseP
#ifdef _CD_TIMING_
Drvg1_CPU = TCpu2-TCpu1
Drvg1_Wall = TWall2-TWall1
write(u6,*) '-------------------------'
write(u6,*) 'Time spent in Prepp:'
write(u6,*) 'Wall/CPU',Prepp_Wall,Prepp_CPU
write(u6,*) '-------------------------'
write(u6,*) 'Time spent in Pget:'
write(u6,*) 'Wall/CPU',Pget_Wall,Pget_CPU
write(u6,*) '-------------------------'
write(u6,*) 'Time spent in Drvg1:'
write(u6,*) 'Wall/CPU',Drvg1_Wall,Drvg1_CPU
write(u6,*) '-------------------------'
Total_Dens_Wall = Prepp_Wall+Pget_Wall
Total_Dens_CPU = Prepp_CPU+Pget_CPU
Total_Der_Wall = Drvg1_Wall-Total_Dens_Wall
Total_Der_CPU = Drvg1_CPU-Total_Dens_CPU
Total_Der_Wall2 = TwoEl_Wall
Total_Der_CPU2 = TwoEl_CPU

write(u6,*) 'Total Time for Density:'
write(u6,*) 'Wall/CPU',Total_Dens_Wall,Total_Dens_CPU
write(u6,*) '-------------------------'
write(u6,*) 'Total Time for Derivatives:'
write(u6,*) 'Wall/CPU',Total_Der_Wall2,Total_Der_CPU2
write(u6,*) '-------------------------'
write(u6,*) 'Derivative check:'
write(u6,*) 'Wall/CPU',Total_Der_Wall,Total_Der_CPU
write(u6,*) '-------------------------'
#endif
Verbose = .false.
FreeK2 = .true.
call Term_Ints(Verbose,FreeK2)
!                                                                      *
!***********************************************************************
!                                                                      *
call Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)

iPren = 3+max(1,int(log10(Pren+1.0e-3_wp)))
iPrem = 3+max(1,int(log10(Prem+1.0e-3_wp)))
write(formt,'(A,I2,A,I2,A)') '(A,F',iPren,'.0,A,F',iPrem,'.0,A)'
if (iPrint >= 6) then
  write(u6,formt) ' A total of',Pren,' entities were prescreened and',Prem,' were kept.'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()

return

end subroutine Drvg1
