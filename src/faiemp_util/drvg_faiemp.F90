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
! Copyright (C) Ben Swerts                                             *
!***********************************************************************

subroutine Drvg_FAIEMP(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
!  Object: driver for the derivatives of central-fragment              *
!          two-electron integrals. The four outermost loops            *
!          will controll the type of the two-electron integral, eg.    *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis func-  *
!          tions of the requested type.                                *
!                                                                      *
!     Author: Ben Swerts                                               *
!                                                                      *
!     based on Drvg1                                                   *
!***********************************************************************

use k2_setup, only: Data_k2
use iSD_data, only: iSD
use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
use Basis_Info, only: nBas, nBas_Frag, Shells
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Eight, Half
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
real(kind=wp) :: Coor(3,4), PMax, A_int, Cnt, P_Eff, Prem, Pren, TCpu1, TCpu2, ThrAO, TMax_all, TskHi, TskLw, TWall1, TWall2
integer(kind=iwp) :: iAnga(4), iCmpa(4), iShela(4), iShlla(4), iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4), iSD4(0:nSD,4), &
                     MemMax, nBas_Valence(0:7), iRout, iPrint, nBT, nBVT, i, j, iAng, iBasi, iBasn, iS, jS, iBasAO, iBsInc, iCar, &
                     ijklA, ijS, iOpt, ijMax, ipEI, ipEta, ipiEta, ipMem1, ipMem2, ipP, ipQ, iPrem, iPren, ipxA, ipxB, ipxG, ipxD, &
                     ipZi, Mem1, Mem2, iPrimi, iPrInc, jAng, iSh, jBasAO, jBasj, jBasn, jBsInc, jPrInc, k2ij, k2kl, jPrimj, &
                     kBasAO, kBasn, kBask, kBsInc, kBtch, klS, kPrimk, kPrInc, kS, lBasAO, lBasl, lBasn, lBsInc, lPriml, lPrInc, &
                     mBtch, lS, mdci, mdcj, mdck, mdcl, MemPSO, nab, ncd, nDCRR, nDCRS, nEta, nHmab, nHmcd, nHrrab, nij, nijkl, &
                     nPairs, nQuad, nRys, nSkal, nSkal_Fragments, nSkal_Valence, nSO, nZeta, nBtch
logical(kind=iwp) :: EQ, Shijij, AeqB, CeqD, lDummy, DoGrad, DoFock, Indexation, JfGrad(3,4), ABCDeq, No_Batch, Rsv_GTList, &
                     FreeK2, Verbose, Triangular, lNoSkip
character(len=72) :: formt
integer(kind=iwp), save :: MemPrm
integer(kind=iwp), allocatable :: ij(:)
real(kind=wp), allocatable :: TMax(:,:)
external :: Rsv_GTList

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 203
iPrint = nPrint(iRout)
iFnc(1) = 0
iFnc(2) = 0
iFnc(3) = 0
iFnc(4) = 0
PMax = Zero

! Handle both the valence and the fragment basis set

call Set_Basis_Mode('Valence')
call Nr_Shells(nSkal_Valence)
call Set_Basis_Mode('WithFragments')
call SetUp_iSD()
nBT = 0
nBVT = 0
do i=0,nIrrep-1
  nBas_Valence(i) = nBas(i)
  nBVT = nBVT+nBas(i)*(nBas(i)+1)/2
  nBas(i) = nBas(i)+nBas_Frag(i)
  nBT = nBT+nBas(i)*(nBas(i)+1)/2
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Precompute k2 entities.

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
nSkal_Fragments = nSkal-nSkal_Valence
if (iPrint >= 99) write(u6,*) 'nSkal, nSkal_Valence, nSkal_Frag = ',nSkal,nSkal_Valence,nSkal_Fragments
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Prepare handling of the combined two-particle density.
!
call PrepP_FAIEMP(nBas_Valence,nBT,nBVT)
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
!---  Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal,label='TMax')
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

call mma_allocate(ij,nSkal*(nSkal+1),label='ij')
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      ij(nij*2-1) = iS
      ij(nij*2) = jS
    end if
  end do
end do
P_Eff = real(nij,kind=wp)
!                                                                      *
!***********************************************************************
!                                                                      *
!-------Compute FLOP's for the transfer equation.

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
call Init_PPList()
call Init_GTList()
iOpt = 0
Temp(:) = Zero
if (iPrint >= 15) call PrGrad(' In Drvg_FAIEMP: Total Grad (1)',Grad,nGrad,ChDisp)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MemMax)
if (MemMax > 1000) MemMax=MemMax-1000
call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
ipMem1 = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
do
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.
  if (.not. Rsv_GTList(TskLw,TskHi,iOpt,lDummy)) exit

  ! Now do a quadruple loop over shells

  ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
  iS = ij(ijS*2-1)
  jS = ij(ijS*2)
  klS = int(TskLw-real(ijS,kind=wp)*(real(ijS,kind=wp)-One)/Two)
  kS = ij(klS*2-1)
  lS = ij(klS*2)
  Cnt = TskLw
  call CWTime(TCpu1,TWall1)
  do
    if (Cnt-TskHi > 1.0e-10_wp) exit

    A_int = TMax(iS,jS)*TMax(kS,lS)

    lNoSkip = A_Int >= CutInt
    ! only calculate needed integrals and only update the valence part of the
    ! Fock matrix (iS > nSkal_Valence, lS <= nSkal_Valence, jS and kS
    ! belonging to different regions)
    if (jS <= nSkal_Valence) then
      lNoSkip = lNoSkip .and. kS > nSkal_Valence
    else
      lNoSkip = lNoSkip .and. kS <= nSkal_Valence
    end if
    lNoSkip = lNoSkip .and. lS <= nSkal_Valence
    if (lNoSkip) then
      !                                                                *
      !*****************************************************************
      !                                                                *
      call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)
      call Size_SO_block_g(iSD4,nSD,nSO,No_batch)
    end if
    if (No_batch) lNoSkip = .false.
    if (lNoSkip) then

      call Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! --------> Memory Managment <--------
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
      if (nIrrep == 1 .and. ABCDeq .and. mod(ijklA,2) == 1 .and. iPrint > 15) write(u6,*) 'ABCDeq & ijklA'
    end if
    if (nIrrep == 1 .and. ABCDeq .and. mod(ijklA,2) == 1) lNoSkip = .false.
    if (lNoSkip) then
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Decide on the partioning of the shells based on the
      ! available memory and the requested memory.
      !
      ! Now check if all blocks can be computed and stored at
      ! once.

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

              !--Get the 2nd order density matrix in SO basis.

              nijkl = iBasn*jBasn*kBasn*lBasn
              call PGet0(iCmpa,iBasn,jBasn,kBasn,lBasn,Shijij,iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,iFnc(1)*iBasn,iFnc(2)*jBasn, &
                         iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
              if (A_Int*PMax >= CutInt) then

                !--Compute gradients of shell quadruplet

                call TwoEl_g(Coor,iAnga,iCmpa,iShela,iShlla,iAOV,mdci,mdcj,mdck,mdcl,nRys,Data_k2(k2ij),nab,nHmab,nDCRR, &
                             Data_k2(k2kl),ncd,nHmcd,nDCRS,Pren,Prem,iPrimi,iPrInc,jPrimj,jPrInc,kPrimk,kPrInc,lPriml,lPrInc, &
                             Shells(iSD4(0,1))%pCff(iPrimi,iBasAO),iBasn,Shells(iSD4(0,2))%pCff(jPrimj,jBasAO),jBasn, &
                             Shells(iSD4(0,3))%pCff(kPrimk,kBasAO),kBasn,Shells(iSD4(0,4))%pCff(lPriml,lBasAO),lBasn, &
                             Mem_DBLE(ipZeta),Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,Mem_DBLE(ipEta),Mem_DBLE(ipEI),Mem_DBLE(ipQ), &
                             nEta,Mem_DBLE(ipxA),Mem_DBLE(ipxB),Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,JfGrad,JndGrd, &
                             Sew_Scr(ipMem1),nSO,Sew_Scr(ipMem2),Mem2,Aux,nAux,Shijij)

                if (iPrint >= 15) call PrGrad(' In Drvg_FAIEMP: Grad',Temp,nGrad,ChDisp)

              end if
            end do
          end do

        end do
      end do

    end if
    Cnt = Cnt+One
    klS = klS+1
    if (klS > ijS) then
      ijS = ijS+1
      klS = 1
    end if
    iS = ij(ijS*2-1)
    jS = ij(ijS*2)
    kS = ij(klS*2-1)
    lS = ij(klS*2)
  end do

  call CWTime(TCpu2,TWall2)
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
call Free_GTList()
call Free_PPList()
call Free_TList()
call mma_deallocate(ij)
call mma_deallocate(TMax)
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

iPren = 3+max(1,int(log10(Pren+0.001_wp)))
iPrem = 3+max(1,int(log10(Prem+0.001_wp)))
write(formt,'(A,I2,A,I2,A)') '(A,F',iPren,'.0,A,F',iPrem,'.0,A)'
if (iPrint >= 6) then
  write(u6,formt) ' A total of',Pren,' entities were prescreened and',Prem,' were kept.'
end if
! Accumulate the final results
call DScal_(nGrad,Half,Temp,1)
if (iPrint >= 15) call PrGrad('The FAIEMP 2-electron Contribution',Temp,nGrad,ChDisp)
call daxpy_(nGrad,One,Temp,1,Grad,1)

call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Drvg_FAIEMP
