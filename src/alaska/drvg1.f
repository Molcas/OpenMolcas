************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1990-1992,2000, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drvg1(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals. The four outermost loops *
*          will control the type of the two-electron integral, eg.     *
*          (ss|ss), (sd|pp), etc. The next four loops will generate    *
*          list of symmetry distinct centers that do have basis        *
*          functions of the requested type.                            *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              SetUp_Ints                                              *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              Swap                                                    *
*              MemRg1                                                  *
*              PSOAO1                                                  *
*              PGet0                                                   *
*              TwoEl                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified for gradient calculation. January '92           *
*             Modified for SetUp_Ints. January '00                     *
************************************************************************
      use k2_setup
      use iSD_data
      use PSO_Stuff
      use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
      use Basis_Info
      use Sizes_of_Seward, only:S
      use Real_Info, only: CutInt
      Implicit Real*8 (A-H,O-Z)
      External Rsv_GTList
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
#include "para_info.fh"
*     Local arrays
      Real*8  Coor(3,4), Grad(nGrad), Temp(nGrad)
      Integer iAnga(4), iCmpa(4), iShela(4),iShlla(4),
     &        iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4)
      Integer nHrrTb(0:iTabMx,0:iTabMx,2)
      Logical EQ, Shijij, AeqB, CeqD, lDummy,
     &        DoGrad, DoFock, Indexation,
     &        JfGrad(3,4), ABCDeq, No_Batch, Rsv_GTList,
     &        FreeK2, Verbose, Triangular
      Character Format*72
      Character*8 Method_chk
      Real*8, Allocatable:: TMax(:,:)
      Integer, Allocatable:: Ind_ij(:,:)
************ columbus interface ****************************************
      Integer  Columbus
*
      Integer iSD4(0:nSD,4)
      save MemPrm
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      iPrint = 000000000

      iFnc(1)=0
      iFnc(2)=0
      iFnc(3)=0
      iFnc(4)=0
      PMax=Zero
      idum=0
      idum1=0
#ifdef _CD_TIMING_
      Twoel_CPU = 0.0d0
      Twoel_Wall = 0.0d0
      Pget_CPU = 0.0d0
      Pget_Wall = 0.0d0
#endif
      Call QEnter('Drvg1')
      call dcopy_(nGrad,[Zero],0,Temp,1)
*
      Call StatusLine(' Alaska:',' Computing 2-electron gradients')
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*-----Precompute k2 entities.
*
      Call Get_iScalar('Columbus',Columbus)
      Indexation=.False.
*     MP2 gradients:
      Call Get_cArray('Relax Method',Method_chk,8)
      If(Method_chk.eq.'MBPT2   ') Indexation=.True.
************ columbus interface ****************************************
* in order to access the half-sorted density matrix file
* some index arrays must be generated
      If (Columbus.eq.1) Indexation=.True.
      DoFock=.False.
      DoGrad=.True.
      ThrAO=Zero
      Call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
      mSkal=nSkal
      nPairs=nSkal*(nSkal+1)/2
      nQuad =nPairs*(nPairs+1)/2
      Pren = Zero
      Prem = Zero
*                                                                      *
************************************************************************
*                                                                      *
*-----Prepare handling of two-particle density.
*
      Call PrepP

*                                                                      *
************************************************************************
*                                                                      *
      MxPrm = 0
      Do iAng = 0, iAngMx
         MxPrm = Max(MxPrm,S%MaxPrm(iAng))
      End Do
      nZeta = MxPrm * MxPrm
      nEta  = MxPrm * MxPrm
*
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
      Call Shell_MxSchwz(nSkal,TMax)
      TMax_all=Zero
      Do iS = 1, nSkal
         Do jS = 1, iS
            TMax_all=Max(TMax_all,TMax(iS,jS))
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of non-vanishing pairs
*
      Call mma_allocate(Ind_ij,2,nskal*(nSkal+1)/2,Label='Ind_ij')
      nij=0
      Do iS = 1, nSkal
         Do jS = 1, iS
            If (TMax_All*TMax(iS,jS).ge.CutInt) Then
               nij = nij + 1
               Ind_ij(1,nij)=iS
               Ind_ij(2,nij)=jS
            End If
         End Do
      End Do
      P_Eff=Dble(nij)
*                                                                      *
************************************************************************
*                                                                      *
*-------Compute FLOPs for the transfer equation.
*
        Do iAng = 0, iAngMx
           Do jAng = 0, iAng
              nHrrab = 0
              Do i = 0, iAng+1
                 Do j = 0, jAng+1
                    If (i+j.le.iAng+jAng+1) Then
                       ijMax = Min(iAng,jAng)+1
                       nHrrab = nHrrab + ijMax*2+1
                    End If
                 End Do
              End Do
              nHrrTb(iAng,jAng,1)=nHrrab
              nHrrTb(jAng,iAng,1)=nHrrab
           End Do
        End Do

*                                                                      *
************************************************************************
*                                                                      *
      Triangular=.True.
      Call Alloc_TList(Triangular,P_Eff)
      Call Init_TList(Triangular,P_Eff)
      Call Init_PPList
      Call Init_GTList
      iOpt=0
*                                                                      *
************************************************************************
*                                                                      *
*     In MPP case dispatch one processor to do 1-el gradients first
*
      If (nProcs.gt.1.and.King()) Then
            Call Drvh1(Grad,Temp,nGrad)
*        If (nPrint(1).ge.15)
*    &   Call PrGrad(' Gradient excluding two-electron contribution',
*    &               Grad,lDisp(0),ChDisp,5)
         call dcopy_(nGrad,[Zero],0,Temp,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_MaxDBLE(MemMax)
      Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
      ipMem1 = 1
*                                                                      *
************************************************************************
*                                                                      *
*     big loop over individual tasks, distributed over individual nodes
   10 Continue
*     make reservation of a task on global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_GTList(TskLw,TskHi,iOpt,lDummy)) Go To 11
*
*     Now do a quadruple loop over shells
*
      ijS = Int((One+sqrt(Eight*TskLw-Three))/Two)
      iS = Ind_ij(1,ijS)
      jS = Ind_ij(2,ijS)
      klS = Int(TskLw-DBLE(ijS)*(DBLE(ijS)-One)/Two)
      kS = Ind_ij(1,klS)
      lS = Ind_ij(2,klS)
      Count=TskLw
      Call CWTime(TCpu1,TWall1)
  13  Continue
*
         Aint=TMax(iS,jS)*TMax(kS,lS)
         If (AInt.lt.CutInt) Go To 14
         If (iPrint.ge.15) Write (6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
*
************************************************************************
*                                                                      *
         Call Gen_iSD4(iS, jS, kS, lS,iSD,nSD,iSD4)
         Call Size_SO_block_g(iSD4,nSD,nSO,No_batch)
         If (No_batch) Go To 140
*
         Call Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)

*
*                                                                      *
************************************************************************
*                                                                      *
*       --------> Memory Managment <--------
*
*        Compute memory request for the primitives, i.e.
*        how much memory is needed up to the transfer
*        equation.
*
         Call MemRys_g(iSD4,nSD,nRys,MemPrm)
*                                                                      *
************************************************************************
*                                                                      *
         ABCDeq=EQ(Coor(1,1),Coor(1,2)) .and.
     &          EQ(Coor(1,1),Coor(1,3)) .and.
     &          EQ(Coor(1,1),Coor(1,4))
         ijklA=iSD4(1,1)+iSD4(1,2)
     &        +iSD4(1,3)+iSD4(1,4)
         If (nIrrep.eq.1.and.ABCDeq.and.Mod(ijklA,2).eq.1)
     &      Go To 140
*                                                                      *
************************************************************************
*                                                                      *
*        Decide on the partioning of the shells based on the
*        available memory and the requested memory.
*
*        Now check if all blocks can be computed and stored at
*        once.
*

         Call SOAO_g(iSD4,nSD,nSO,
     &               MemPrm, MemMax,
     &               iBsInc,jBsInc,kBsInc,lBsInc,
     &               iPrInc,jPrInc,kPrInc,lPrInc,
     &               ipMem1,ipMem2, Mem1,  Mem2,
     &               iPrint,iFnc, MemPSO)
         iBasi    = iSD4(3,1)
         jBasj    = iSD4(3,2)
         kBask    = iSD4(3,3)
         lBasl    = iSD4(3,4)
*                                                                      *
************************************************************************
*                                                                      *
         Call Int_Parm_g(iSD4,nSD,iAnga,
     &                 iCmpa,iShlla,iShela,
     &                 iPrimi,jPrimj,kPrimk,lPriml,
     &                 indij,k2ij,nDCRR,k2kl,nDCRS,
     &                 mdci,mdcj,mdck,mdcl,AeqB,CeqD,
     &                 nZeta,nEta,ipZeta,ipZI,
     &                 ipP,ipEta,ipEI,ipQ,ipiZet,ipiEta,
     &                 ipxA,ipxB,ipxG,ipxD,l2DI,nab,nHmab,ncd,nHmcd,
     &                 nIrrep)
*                                                                      *
************************************************************************
*                                                                      *
*        Scramble arrays (follow angular index)
*

         Do iCar = 1, 3
            Do iSh = 1, 4
               JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
               If (iAnd(iSD4(15,iSh),2**(iCar-1)) .eq.
     &             2**(iCar-1)) Then
                  JfGrad(iCar,iSh) = .True.
               Else
                  JfGrad(iCar,iSh) = .False.
               End If
            End Do
         End Do
*
         Do 400 iBasAO = 1, iBasi, iBsInc
           iBasn=Min(iBsInc,iBasi-iBasAO+1)
           iAOst(1) = iBasAO-1
         Do 410 jBasAO = 1, jBasj, jBsInc
           jBasn=Min(jBsInc,jBasj-jBasAO+1)
           iAOst(2) = jBasAO-1
         Do 420 kBasAO = 1, kBask, kBsInc
           kBasn=Min(kBsInc,kBask-kBasAO+1)
           iAOst(3) = kBasAO-1
         Do 430 lBasAO = 1, lBasl, lBsInc
           lBasn=Min(lBsInc,lBasl-lBasAO+1)
           iAOst(4) = lBasAO-1
*
*----------Get the 2nd order density matrix in SO basis.
*
           nijkl = iBasn*jBasn*kBasn*lBasn


*     Fetch the T_i,j,kappa, lambda corresponding to
*     kappa = k, lambda = l
*
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU1,Pget0WALL1)
#endif
           Call PGet0(iCmpa,
     &                iBasn,jBasn,kBasn,lBasn,Shijij,
     &                iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,
     &                iFnc(1)*iBasn,iFnc(2)*jBasn,
     &                iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,
     &                Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
           If (AInt*PMax.lt.CutInt) Go To 430
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU2,Pget0WALL2)
           Pget_CPU = Pget_CPU + Pget0CPU2-Pget0CPU1
           Pget_Wall = Pget_Wall + Pget0WALL2-Pget0WALL1
#endif
           If (AInt*PMax.lt.CutInt) Go To 430
*
*----------Compute gradients of shell quadruplet
*
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU1,TwoelWall1) ! timing_cdscf
#endif
           Call TwoEl_g(Coor,
     &          iAnga,iCmpa,iShela,iShlla,iAOV,
     &          mdci,mdcj,mdck,mdcl,nRys,
     &          Data_k2(k2ij),nab,nHmab,nDCRR,
     &          Data_k2(k2kl),ncd,nHmcd,nDCRS,Pren,Prem,
     &          iPrimi,iPrInc,jPrimj,jPrInc,
     &          kPrimk,kPrInc,lPriml,lPrInc,
     &          Shells(iSD4(0,1))%pCff(1,iBasAO),iBasn,
     &          Shells(iSD4(0,2))%pCff(1,jBasAO),jBasn,
     &          Shells(iSD4(0,3))%pCff(1,kBasAO),kBasn,
     &          Shells(iSD4(0,4))%pCff(1,lBasAO),lBasn,
     &          Mem_DBLE(ipZeta),Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,
     &          Mem_DBLE(ipEta), Mem_DBLE(ipEI),Mem_DBLE(ipQ),nEta,
     &          Mem_DBLE(ipxA),Mem_DBLE(ipxB),
     &          Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,
     &          JfGrad,JndGrd,Sew_Scr(ipMem1), nSO,Sew_Scr(ipMem2),Mem2,
     &          Aux,nAux,Shijij)
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU2,TwoelWall2)
           Twoel_CPU = Twoel_CPU + TwoelCPU2-TwoelCPU1
           Twoel_Wall = Twoel_Wall + TwoelWall2-TwoelWall1
#endif
            If (iPrint.ge.15)
     &         Call PrGrad(' In Drvg1: Grad',Temp,nGrad,ChDisp,5)
*
 430     Continue
 420     Continue
*
 410     Continue
 400     Continue
*
 140     Continue
*
 14      Continue
         Count=Count+One
         If (Count-TskHi.gt.1.0D-10) Go To 12
         klS = klS + 1
         If (klS.gt.ijS) Then
            ijS = ijS + 1
            klS = 1
         End If
         iS = Ind_ij(1,ijS)
         jS = Ind_ij(2,ijS)
         kS = Ind_ij(1,klS)
         lS = Ind_ij(2,klS)
         Go To 13
*
*     Task endpoint
 12   Continue
      Call CWTime(TCpu2,TWall2)
      Call SavTim(4,TCpu2-TCpu1,TWall2-Twall1)
      Call SavStat(1,One,'+')
      Call SavStat(2,TskHi-TskLw+One,'+')
      Go To 10
 11   Continue
*     End of big task loop
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _MOLCAS_MPP_
      Call GADGOP(Temp,nGrad,'+')
#endif
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Sew_Scr)
      Call Free_GTList
      Call Free_PPList
      Call Free_TList
      Call mma_deallocate(Ind_ij)
      Call mma_deallocate(TMax)
*                                                                      *
************************************************************************
*                                                                      *
      Call CloseP
#ifdef _CD_TIMING_
      Drvg1_CPU = TCpu2-TCpu1
      Drvg1_Wall= TWall2-TWall1
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Prepp:'
      Write(6,*) 'Wall/CPU',Prepp_Wall, Prepp_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Pget:'
      Write(6,*) 'Wall/CPU',Pget_Wall, Pget_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Drvg1:'
      Write(6,*) 'Wall/CPU',Drvg1_Wall, Drvg1_CPU
      Write(6,*) '-------------------------'
      Total_Dens_Wall = Prepp_Wall+Pget_Wall
      Total_Dens_CPU = Prepp_CPU+Pget_CPU
      Total_Der_Wall = Drvg1_Wall - Total_Dens_Wall
      Total_Der_CPU = Drvg1_CPU - Total_Dens_CPU
      Total_Der_Wall2 = TwoEl_Wall
      Total_Der_CPU2 = TwoEl_CPU

      Write(6,*) 'Total Time for Density:'
      Write(6,*) 'Wall/CPU',Total_Dens_Wall, Total_Dens_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Total Time for Derivatives:'
      Write(6,*) 'Wall/CPU',Total_Der_Wall2, Total_Der_CPU2
      Write(6,*) '-------------------------'
      Write(6,*) 'Derivative check:'
      Write(6,*) 'Wall/CPU',Total_Der_Wall, Total_Der_CPU
      Write(6,*) '-------------------------'
#endif
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*                                                                      *
************************************************************************
*                                                                      *
*
      Call Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)
*
      iPren=3+Max(1,Int(Log10(Pren+0.001D+00)))
      iPrem=3+Max(1,Int(Log10(Prem+0.001D+00)))
      Write (Format,'(A,I2,A,I2,A)') '(A,F',iPren,
     &           '.0,A,F',iPrem,'.0,A)'
      If (iPrint.ge.6) Then
      Write (6,Format)
     &   ' A total of', Pren,' entities were prescreened and',
     &                  Prem,' were kept.'
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      Call QExit('Drvg1')
      Return
      End
