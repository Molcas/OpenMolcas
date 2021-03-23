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

use k2_setup
use iSD_data
use PSO_Stuff
use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
use Basis_Info
use Sizes_of_Seward, only: S
use Real_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Para_Info, only: nProcs, King

implicit real*8(A-H,O-Z)
external Rsv_GTList
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
! Local arrays
real*8 Coor(3,4), Grad(nGrad), Temp(nGrad)
integer iAnga(4), iCmpa(4), iShela(4), iShlla(4), iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4)
logical EQ, Shijij, AeqB, CeqD, lDummy, DoGrad, DoFock, Indexation, JfGrad(3,4), ABCDeq, No_Batch, Rsv_GTList, FreeK2, Verbose, &
  Triangular
character format*72
character*8 Method_chk
real*8, allocatable :: TMax(:,:)
integer, allocatable :: Ind_ij(:,:)
!*********** columbus interface ****************************************
integer Columbus

integer iSD4(0:nSD,4)
save MemPrm

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 9
iPrint = nPrint(iRout)
iPrint = 000000000

iFnc(1) = 0
iFnc(2) = 0
iFnc(3) = 0
iFnc(4) = 0
PMax = Zero
#ifdef _CD_TIMING_
Twoel_CPU = 0.0d0
Twoel_Wall = 0.0d0
Pget_CPU = 0.0d0
Pget_Wall = 0.0d0
#endif
call dcopy_(nGrad,[Zero],0,Temp,1)

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
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      Ind_ij(1,nij) = iS
      Ind_ij(2,nij) = jS
    end if
  end do
end do
P_Eff = dble(nij)
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
  !if (nPrint(1) >= 15) call PrGrad(' Gradient excluding two-electron contribution',Grad,lDisp(0),ChDisp,5)
  call dcopy_(nGrad,[Zero],0,Temp,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MemMax)
call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
ipMem1 = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
10 continue
! make reservation of a task on global task list and get task range
! in return. Function will be false if no more tasks to execute.
if (.not. Rsv_GTList(TskLw,TskHi,iOpt,lDummy)) Go To 11

! Now do a quadruple loop over shells

ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
iS = Ind_ij(1,ijS)
jS = Ind_ij(2,ijS)
klS = int(TskLw-dble(ijS)*(dble(ijS)-One)/Two)
kS = Ind_ij(1,klS)
lS = Ind_ij(2,klS)
Count = TskLw
call CWTime(TCpu1,TWall1)
13 continue

Aint = TMax(iS,jS)*TMax(kS,lS)
if (AInt < CutInt) Go To 14
if (iPrint >= 15) write(6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
!                                                                      *
!***********************************************************************
!                                                                      *
call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)
call Size_SO_block_g(iSD4,nSD,nSO,No_batch)
if (No_batch) Go To 140

call Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)
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
if ((nIrrep == 1) .and. ABCDeq .and. (mod(ijklA,2) == 1)) Go To 140
!                                                                      *
!***********************************************************************
!                                                                      *
! Decide on the partioning of the shells based on the
! available memory and the requested memory.
!
! Now check if all blocks can be computed and stored at once.

call SOAO_g(iSD4,nSD,nSO,MemPrm,MemMax,iBsInc,jBsInc,kBsInc,lBsInc,iPrInc,jPrInc,kPrInc,lPrInc,ipMem1,ipMem2,Mem1,Mem2,iPrint, &
            iFnc,MemPSO)
iBasi = iSD4(3,1)
jBasj = iSD4(3,2)
kBask = iSD4(3,3)
lBasl = iSD4(3,4)
!                                                                      *
!***********************************************************************
!                                                                      *
call Int_Parm_g(iSD4,nSD,iAnga,iCmpa,iShlla,iShela,iPrimi,jPrimj,kPrimk,lPriml,indij,k2ij,nDCRR,k2kl,nDCRS,mdci,mdcj,mdck,mdcl, &
                AeqB,CeqD,nZeta,nEta,ipZeta,ipZI,ipP,ipEta,ipEI,ipQ,ipiZet,ipiEta,ipxA,ipxB,ipxG,ipxD,l2DI,nab,nHmab,ncd,nHmcd, &
                nIrrep)
!                                                                      *
!***********************************************************************
!                                                                      *
! Scramble arrays (follow angular index)

do iCar=1,3
  do iSh=1,4
    JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
    if (iand(iSD4(15,iSh),2**(iCar-1)) == 2**(iCar-1)) then
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

#       ifdef _CD_TIMING_
        call CWTIME(Pget0CPU1,Pget0WALL1)
#       endif
        call PGet0(iCmpa,iBasn,jBasn,kBasn,lBasn,Shijij,iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,iFnc(1)*iBasn,iFnc(2)*jBasn, &
                   iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
        if (AInt*PMax < CutInt) Go To 430
#       ifdef _CD_TIMING_
        call CWTIME(Pget0CPU2,Pget0WALL2)
        Pget_CPU = Pget_CPU+Pget0CPU2-Pget0CPU1
        Pget_Wall = Pget_Wall+Pget0WALL2-Pget0WALL1
#       endif
        if (AInt*PMax < CutInt) Go To 430

        ! Compute gradients of shell quadruplet

#       ifdef _CD_TIMING_
        call CWTIME(TwoelCPU1,TwoelWall1) ! timing_cdscf
#       endif
        call TwoEl_g(Coor,iAnga,iCmpa,iShela,iShlla,iAOV,mdci,mdcj,mdck,mdcl,nRys,Data_k2(k2ij),nab,nHmab,nDCRR,Data_k2(k2kl),ncd, &
                     nHmcd,nDCRS,Pren,Prem,iPrimi,iPrInc,jPrimj,jPrInc,kPrimk,kPrInc,lPriml,lPrInc, &
                     Shells(iSD4(0,1))%pCff(1,iBasAO),iBasn,Shells(iSD4(0,2))%pCff(1,jBasAO),jBasn, &
                     Shells(iSD4(0,3))%pCff(1,kBasAO),kBasn,Shells(iSD4(0,4))%pCff(1,lBasAO),lBasn,Mem_DBLE(ipZeta), &
                     Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,Mem_DBLE(ipEta),Mem_DBLE(ipEI),Mem_DBLE(ipQ),nEta,Mem_DBLE(ipxA), &
                     Mem_DBLE(ipxB),Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,JfGrad,JndGrd,Sew_Scr(ipMem1),nSO,Sew_Scr(ipMem2), &
                     Mem2,Aux,nAux,Shijij)
#       ifdef _CD_TIMING_
        call CWTIME(TwoelCPU2,TwoelWall2)
        Twoel_CPU = Twoel_CPU+TwoelCPU2-TwoelCPU1
        Twoel_Wall = Twoel_Wall+TwoelWall2-TwoelWall1
#       endif
        if (iPrint >= 15) call PrGrad(' In Drvg1: Grad',Temp,nGrad,ChDisp,5)

430     continue
      end do
    end do

  end do
end do

140 continue

14 continue
Count = Count+One
if (Count-TskHi > 1.0D-10) Go To 12
klS = klS+1
if (klS > ijS) then
  ijS = ijS+1
  klS = 1
end if
iS = Ind_ij(1,ijS)
jS = Ind_ij(2,ijS)
kS = Ind_ij(1,klS)
lS = Ind_ij(2,klS)
Go To 13

! Task endpoint
12 continue
call CWTime(TCpu2,TWall2)
call SavTim(4,TCpu2-TCpu1,TWall2-Twall1)
call SavStat(1,One,'+')
call SavStat(2,TskHi-TskLw+One,'+')
Go To 10
11 continue
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
!                                                                      *
!***********************************************************************
!                                                                      *
call CloseP
#ifdef _CD_TIMING_
Drvg1_CPU = TCpu2-TCpu1
Drvg1_Wall = TWall2-TWall1
write(6,*) '-------------------------'
write(6,*) 'Time spent in Prepp:'
write(6,*) 'Wall/CPU',Prepp_Wall,Prepp_CPU
write(6,*) '-------------------------'
write(6,*) 'Time spent in Pget:'
write(6,*) 'Wall/CPU',Pget_Wall,Pget_CPU
write(6,*) '-------------------------'
write(6,*) 'Time spent in Drvg1:'
write(6,*) 'Wall/CPU',Drvg1_Wall,Drvg1_CPU
write(6,*) '-------------------------'
Total_Dens_Wall = Prepp_Wall+Pget_Wall
Total_Dens_CPU = Prepp_CPU+Pget_CPU
Total_Der_Wall = Drvg1_Wall-Total_Dens_Wall
Total_Der_CPU = Drvg1_CPU-Total_Dens_CPU
Total_Der_Wall2 = TwoEl_Wall
Total_Der_CPU2 = TwoEl_CPU

write(6,*) 'Total Time for Density:'
write(6,*) 'Wall/CPU',Total_Dens_Wall,Total_Dens_CPU
write(6,*) '-------------------------'
write(6,*) 'Total Time for Derivatives:'
write(6,*) 'Wall/CPU',Total_Der_Wall2,Total_Der_CPU2
write(6,*) '-------------------------'
write(6,*) 'Derivative check:'
write(6,*) 'Wall/CPU',Total_Der_Wall,Total_Der_CPU
write(6,*) '-------------------------'
#endif
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
call Free_iSD()

return

end subroutine Drvg1
