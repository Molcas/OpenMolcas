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
* Copyright (C) 1990,1991,1992,2000,2007, Roland Lindh                 *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drvg1_2Center_RI(Grad,Temp,nGrad,ij2,nij_Eff)
************************************************************************
*                                                                      *
*  Object: driver for 2-center two-electron integrals in the RI scheme.*
*                                                                      *
*   The integral derivative is formulated as                           *
*   -Sum(ML) X_ij^K   V_LM^(1) X_kl^L  where                           *
*                                                                      *
*  X_ij^K = Sum(L) R_ij_L  Q_L^K                                       *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified for gradient calculation. January '92           *
*             Modified for SetUp_Ints. January '00                     *
*             Modified for 2-center RI gradients, January '07          *
************************************************************************
      use k2_setup
      use iSD_data
      use pso_stuff
      use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
      use Basis_Info
      use Sizes_of_Seward, only:S
      use Gateway_Info, only: CutInt
      use RICD_Info, only: Do_RI
      use Symmetry_Info, only: nIrrep
      use Para_Info, only: nProcs, King
      use ExTerm, only: CijK, AMP2, iMP2prpt, A
      Implicit Real*8 (A-H,O-Z)
      External Rsv_Tsk
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
#include "exterm.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
      Integer nGrad, nij_Eff
      Real*8  Grad(nGrad), Temp(nGrad)
      Integer, Allocatable :: ij2(:,:)

*     Local arrays
      Real*8  Coor(3,4)
      Integer iAnga(4), iCmpa(4), iShela(4),iShlla(4),
     &        iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4)
      Logical EQ, Shijij, AeqB, CeqD,
     &        DoGrad, DoFock, Indexation, FreeK2, Verbose,
     &        JfGrad(3,4), ABCDeq, No_Batch, Rsv_Tsk
      Character Format*72
*
      Character*8 Method_chk
      ! Character*4096 RealName
*
      Integer iSD4(0:nSD,4)
      Save MemPrm

      Real*8, Allocatable:: TMax2(:,:), TMax1(:), Tmp(:,:)
      Integer, Allocatable:: Shij(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
#ifdef _CD_TIMING_
      Twoel2_CPU = 0.0d0
      Twoel2_Wall = 0.0d0
      Pget2_CPU = 0.0d0
      Pget2_Wall = 0.0d0
#endif

      iFnc(1)=0
      iFnc(2)=0
      iFnc(3)=0
      iFnc(4)=0
      PMax=Zero
      Temp(:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     Handle only the auxiliary basis set.
*
      If (Do_RI) Then
         Call Set_Basis_Mode('Auxiliary')
      Else
         Call Set_Basis_Mode('Valence')
      End If
      Call Setup_iSD()

************************************************************************
*                                                                      *
*-----Precompute k2 entities.
*
      Indexation=.False.
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
      MxPrm = 0
      Do iAng = 0, S%iAngMx
         MxPrm = Max(MxPrm,S%MaxPrm(iAng))
      End Do
      nZeta = MxPrm * MxPrm
      nEta  = MxPrm * MxPrm
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      If (Do_RI) Then
         nTMax=nSkal
         Call mma_allocate(TMax1,nTMax,Label='TMax1')
         Call mma_allocate(Tmp,nSkal,nSkal,Label='Tmp')
         Call Shell_MxSchwz(nSkal,Tmp)
         TMax1(1:nSkal)=Tmp(1:nSkal,nSkal)
         Call mma_deallocate(Tmp)

         TMax_all=Zero
         Do iS = 1, nSkal-1
            TMax_all=Max(TMax_all,TMax1(iS))
         End Do
      Else
         Call mma_allocate(TMax2,nSkal,nSkal,Label='TMax2')
         Call Shell_MxSchwz(nSkal,TMax2)
         TMax_all=Zero
         Do ij = 1, nij_Eff
            iS = ij2(1,ij)
            jS = ij2(2,ij)
            TMax_all=Max(TMax_all,TMax2(iS,jS))
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate some scratch arrays to be used by the pget routines.
*     In particular we will have temporary arrays for A_IJ and C_ijK.
*
*     Lower case: valence basis set
*     Upper case: auxiliary basis sets
*
      If(DoCholExch) Then
*
*        Find the largest number of contractions in any given shell
*        of auxiliary functions.

         MxChVInShl = 1
         If(Do_RI) Then
            Do i = 1, nSkal
               MxChVInShl = max(MxChVInShl,iSD(3,i))
            End Do
         Else
            Write (6,*) 'Not Implemented for Cholesky yet!'
            Call Abend()
         End If
*
*        Scratch for A_IJ
*
         lA = MxChVInShl*MxChVInShl
         Call mma_allocate(A,lA,Label='A')
         If (iMP2Prpt.eq.2) Then
            lA_MP2=MxChVInShl
            Call mma_allocate(AMP2,lA_MP2,2,Label='AMP2')
         End If
*
*        Find the largest set of ij. The basis i and j is due to the
*        CD of the one-particle density.
*
         nIJRMax = 0
         Do jDen = 1,nKvec
            Do iSym1 = 1, nIrrep
               Do iSym2 = 1, nIrrep
                  nIJRMax = max(nIJRMax,nIJR(iSym1,iSym2,jDen))
               End Do
            End Do
         End Do
*
*        Get scratch for C_kl^I and C_kl^J.
*        Note that we need nDen arrays for C_kl^I and one for C_kl^J
*        A_IJ = Sum(kl) C_kl^I x C_kl^J
*
         Call mma_allocate(CijK,nIJRMax*MxChVInShl*(nKvec+1),
     &                     Label='CijK')
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of non-vanishing pairs
*
      If (Do_RI) Then
         mij=(nSkal-1)
         Call mma_allocate(Shij,2,mij,Label='Shij')
         nij=0
         Do iS = 1, nSkal-1
            If (TMax_All*TMax1(iS).ge.CutInt) Then
               nij = nij + 1
               Shij(1,nij)=nSkal
               Shij(2,nij)=iS
            End If
         End Do
      Else
         mij=nij_Eff
         Call mma_allocate(Shij,2,mij,Label='Shij')
         nij=0
         Do ij = 1, nij_Eff
            iS = ij2(1,ij)
            jS = ij2(2,ij)
            If (TMax_All*TMax2(iS,jS).ge.CutInt) Then
               nij = nij + 1
               Shij(1,nij)=iS
               Shij(2,nij)=jS
            End If
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----CASPT2
*
      Call Get_cArray('Relax Method',Method_chk,8)
      If (Method_chk.eq.'CASPT2  ') Then
        !! Just read A_{JK} type matrix constructed in CASPT2
        nBasT   = 0
        nBasA   = 0
        nBasASQ = 0
        Do iSym = 0, nIrrep-1
          nBasT   = nBasT   + nBas(iSym)
          nBasA   = nBasA   + nBas_Aux(iSym)-1
          nBasASQ = nBasASQ + (nBas_Aux(iSym)-1)**2
        End Do
        Call MMA_Allocate(A_PT2,nBasA,nBasA,Label='A_PT2')
        !! Now, read
    !     Call PrgmTranslate('CMOPT2',RealName,lRealName)
    !     LuCMOPT2 = 61
    !     Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
    !  &                        'DIRECT','UNFORMATTED',
    !  &                        iost,.FALSE.,
    !  &                        1,'OLD',is_error)
    !     Read (LuCMOPT2) A_PT2(1:nBasASQ,1)
    !     Close (LuCMOPT2)
        ! Read A_PT2 from LUAPT2
        LUAPT2 = 77
        call daname_mf_wa(LUAPT2, 'A_PT2')
        id = 0
        call ddafile(LUAPT2, 2, A_PT2, nBasASq, id)
        call daclos(LUAPT2)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute FLOP's for the transfer equation.
*
      Do iAng = 0, S%iAngMx
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
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     For a parallel implementation the iterations over shell-pairs
*     are parallelized.

      Call Init_Tsk(id,nij*(nij+1)/2)
*                                                                      *
************************************************************************
*                                                                      *
*     In MPP case dispatch one processor to do 1-el gradients first
*
      If (nProcs.gt.1.and.King()) Then
         If (Do_RI) Call Free_iSD()
         Call Drvh1(Grad,Temp,nGrad)
*        If (nPrint(1).ge.15)
*    &   Call PrGrad(' Gradient excluding two-electron contribution',
*    &               Grad,lDisp(0),ChDisp)
         call dcopy_(nGrad,[Zero],0,Temp,1)
         If (Do_RI) Then
            Call Set_Basis_Mode('Auxiliary')
            Call Setup_iSD()
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_MaxDBLE(MemMax)
      If (MemMax.gt.1000) MemMax=MemMax-1000
      Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
      ipMem1=1
*                                                                      *
************************************************************************
*                                                                      *
C     If (MyRank.ne.0) Go To 11
*     big loop over individual tasks, distributed over individual nodes
   10 Continue
*     make reservation of a task on global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_Tsk(id,jlS)) Go To 11
*
*     Now do a quadruple loop over shells
*
      jS_= Int((One+sqrt(Eight*DBLE(jlS)-Three))/Two)
      iS = Shij(1,jS_)
      jS = Shij(2,jS_)
      lS_= Int(DBLE(jlS)-DBLE(jS_)*(DBLE(jS_)-One)/Two)
      kS = Shij(1,lS_)
      lS = Shij(2,lS_)
      Call CWTime(TCpu1,TWall1)
*
         If (Do_RI) Then
            Aint=TMax1(jS)*TMax1(lS)
         Else
            Aint=TMax2(iS,jS)*TMax2(kS,lS)
         End If
         If (AInt.lt.CutInt) Go To 10
C        If (is.eq.3.and.js.eq.3.and.ks.eq.1.and.ls.eq.1) Then
C           iPrint=15
C           nPrint(39)=15
C        Else
C           iPrint=nPrint(iRout)
C           nPrint(39)=5
C        End If
         If (iPrint.ge.15) Write (6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
*                                                                      *
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
     &               iFnc, MemPSO)
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
     &                 k2ij,nDCRR,k2kl,nDCRS,
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
               If ((iSh.eq.1 .or. iSh.eq.3).and.Do_RI) Then
                  JfGrad(iCar,iSh) = .False.
               Else If (iAnd(iSD4(15,iSh),2**(iCar-1)) .eq.
     &                                    2**(iCar-1)) Then
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
*
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
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU1,Pget0WALL1)
#endif
           Call PGet0(iCmpa,
     &                iBasn,jBasn,kBasn,lBasn,Shijij,
     &                iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,
     &                iFnc(1)*iBasn,iFnc(2)*jBasn,
     &                iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,
     &                Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU2,Pget0WALL2)
           Pget2_CPU = Pget2_CPU + Pget0CPU2-Pget0CPU1
           Pget2_Wall = Pget2_Wall + Pget0WALL2-Pget0WALL1
#endif
            If (AInt*PMax.lt.CutInt) Go To 430
*
*----------Compute gradients of shell quadruplet
*
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU1,TwoelWall1)
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
           Twoel2_CPU = Twoel2_CPU + TwoelCPU2-TwoelCPU1
           Twoel2_Wall = Twoel2_Wall + TwoelWall2-TwoelWall1
#endif
            If (iPrint.ge.15)
     &         Call PrGrad(' In Drvg1_2Center_RI: Grad',
     &                  Temp,nGrad,ChDisp)
*
 430     Continue
 420     Continue
*
 410     Continue
 400     Continue
*
 140     Continue
*
         Go To 10
 11   Continue
*     End of big task loop
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Sew_Scr)
      Call Free_Tsk(id)
      If (Method_chk.eq.'CASPT2  ') Call MMA_DeAllocate(A_PT2)
      Call mma_deallocate(Shij)
      If (Allocated(TMax1)) Call mma_deallocate(TMax1)
      If (Allocated(TMax2)) Call mma_deallocate(TMax2)
*                                                                      *
************************************************************************
*                                                                      *
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*                                                                      *
************************************************************************
*                                                                      *
      Call Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)
*
      iPren=3+Max(1,Int(Log10(Pren+0.001D+00)))
      iPrem=3+Max(1,Int(Log10(Prem+0.001D+00)))
      Write (Format,'(A,I2,A,I2,A)') '(A,F',iPren,
     &           '.0,A,F',iPrem,'.0,A)'
      If (iPrint.ge.6) Then
         Write (6,Format)
     &      ' A total of', Pren,' entities were prescreened and',
     &                     Prem,' were kept.'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If(DoCholExch) Then
         Call mma_deallocate(CijK)
         Call mma_deallocate(A)
      End If
      If (Allocated(AMP2)) Call mma_deallocate(AMP2)
*
      Call Free_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
