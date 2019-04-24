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
      SubRoutine Drvg1_2Center_RI(Grad,Temp,nGrad,ip_ij2,nij_Eff)
************************************************************************
*                                                                      *
*  Object: driver for 2-center two-electron integrals in the RI scheme.*
*                                                                      *
*   The integral derivative is formulated as                           *
*   -Sum(ML) X_ij^K   V_LM^(1) X_kl^L  where                           *
*                                                                      *
*  X_ij^K = Sum(L) R_ij_L  Q_L^K                                       *
*                                                                      *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              SetUp_Ints                                              *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
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
*             Modified for 2-center RI gradients, January '07          *
*                                                                      *
************************************************************************
      use k2_setup
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
      External Rsv_Tsk
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
#include "exterm.fh"
#include "chomp2g_alaska.fh"
#include "pso.fh"
#include "para_info.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
*     Local arrays
      Real*8  Coor(3,4), Grad(nGrad), Temp(nGrad)
      Integer iAnga(4), iCmpa(4), iShela(4),iShlla(4),
     &        iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4)
      Integer nHrrTb(0:iTabMx,0:iTabMx,2)
      Logical EQ, Shijij, AeqB, CeqD,
     &        DoGrad, DoFock, Indexation, FreeK2, Verbose,
     &        JfGrad(3,4), ABCDeq, No_Batch, Rsv_Tsk
      Character Format*72
*
      Integer iSD4(0:nSD,4)
      Save MemPrm
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      TMax1(i)=Work(ipTMax-1+i)
      TMax2(i,j)=Work(ipTMax-1+(j-1)*nSkal+i)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      Call QEnter('Drvg1_2Center_RI')
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
      idum=0
      idum1=0
      call dcopy_(nGrad,[Zero],0,Temp,1)
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
      Do iAng = 0, iAngMx
         MxPrm = Max(MxPrm,MaxPrm(iAng))
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
         Call GetMem('TMax','Allo','Real',ipTMax,nTMax)
         Call Allocate_Work(ip_Tmp,nSkal**2)
         Call Shell_MxSchwz(nSkal,Work(ip_Tmp))
        call dcopy_(nSkal,Work(ip_Tmp+(nSkal-1)*nSkal),1,Work(ipTMax),1)
         Call Free_Work(ip_Tmp)
         TMax_all=Zero
         Do iS = 1, nSkal-1
            TMax_all=Max(TMax_all,TMax1(iS))
         End Do
      Else
         nTMax=nSkal**2
         Call GetMem('TMax','Allo','Real',ipTMax,nTMax)
         Call Shell_MxSchwz(nSkal,Work(ipTMax))
         TMax_all=Zero
         Do ij = 1, nij_Eff
            iS = iWork(ip_ij2-1 +(ij-1)*2 + 1)
            jS = iWork(ip_ij2-1 +(ij-1)*2 + 2)
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
         Call GetMem('A','Allo','Real',ip_A,lA)
         If (iMP2Prpt.eq.2) Then
            lA_MP2=MxChVInShl
            Call GetMem('A_MP2(1)','Allo','Real',ip_A_MP2(1),lA_MP2)
            Call GetMem('A_MP2(2)','Allo','Real',ip_A_MP2(2),lA_MP2)
         Else
            ip_A_MP2(1) = ip_Dummy
            ip_A_MP2(2) = ip_Dummy
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
         lCijK = nIJRMax*MxChVInShl
         Call GetMem('CijK','Allo','Real',ip_CijK,(nKvec+1)*lCijK)
*
      Else
*
         lCijK = 0
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of non-vanishing pairs
*
      If (Do_RI) Then
         mij=(nSkal-1)*2
         Call GetMem('ip_ij','Allo','Inte',ip_ij,mij)
         nij=0
         Do iS = 1, nSkal-1
            If (TMax_All*TMax1(iS).ge.CutInt) Then
               nij = nij + 1
               iWork(ip_ij + 2*(nij-1)  )=nSkal
               iWork(ip_ij + 2*(nij-1)+1)=iS
            End If
         End Do
      Else
         mij=2*nij_Eff
         Call GetMem('ip_ij','Allo','Inte',ip_ij,mij)
         nij=0
         Do ij = 1, nij_Eff
            iS = iWork(ip_ij2-1 +(ij-1)*2 + 1)
            jS = iWork(ip_ij2-1 +(ij-1)*2 + 2)
            If (TMax_All*TMax2(iS,jS).ge.CutInt) Then
               nij = nij + 1
               iWork((nij-1)*2+ip_ij  )=iS
               iWork((nij-1)*2+ip_ij+1)=jS
            End If
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute FLOP's for the transfer equation.
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
*    &               Grad,lDisp(0),lIrrep,ChDisp,iPrint)
         call dcopy_(nGrad,[Zero],0,Temp,1)
         If (Do_RI) Then
            Call Set_Basis_Mode('Auxiliary')
            Call Setup_iSD()
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('MemMax','Max','Real',iDum,MemMax)
      Call GetMem('MemMax','Allo','Real',ipMem1,MemMax)
*                                                                      *
************************************************************************
*                                                                      *
C     Call Get_MyRank(MyRank_)
C     If (MyRank_.ne.0) Go To 11
*     big loop over individual tasks, distributed over individual nodes
   10 Continue
*     make reservation of a task on global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_Tsk(id,jlS)) Go To 11
*
*     Now do a quadruple loop over shells
*
      jS_= Int((One+sqrt(Eight*DBLE(jlS)-Three))/Two)
      iS = iWork((jS_-1)*2  +ip_ij)
      jS = iWork((jS_-1)*2+1+ip_ij)
      lS_= Int(DBLE(jlS)-DBLE(jS_)*(DBLE(jS_)-One)/Two)
      kS = iWork((lS_-1)*2  +ip_ij)
      lS = iWork((lS_-1)*2+1+ip_ij)
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
         Call Size_SO_block_g(iSD4,nSD,Petite,nSO,No_batch)
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
     &               nExp,nBasis,MxShll,
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
     &                 ipCffi,jpCffj,kpCffk,lpCffl,
     &                 nExp,ipExp,ipCff,MxShll,
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
           Call PGet0(iCmpa,iShela,
     &                iBasn,jBasn,kBasn,lBasn,Shijij,
     &                iAOV,iAOst,nijkl,Work(ipMem1),nSO,
     &                iFnc(1)*iBasn,iFnc(2)*jBasn,
     &                iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,
     &                ipMem2,iS,jS,kS,lS,nQuad,PMax)
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
     &          Work(ipCffi+(iBasAO-1)*iPrimi),iBasn,
     &          Work(jpCffj+(jBasAO-1)*jPrimj),jBasn,
     &          Work(kpCffk+(kBasAO-1)*kPrimk),kBasn,
     &          Work(lpCffl+(lBasAO-1)*lPriml),lBasn,
     &          Work(ipZeta),Work(ipZI),Work(ipP),nZeta,
     &          Work(ipEta), Work(ipEI),Work(ipQ),nEta,
     &          Work(ipxA),Work(ipxB),Work(ipxG),Work(ipxD),Temp,nGrad,
     &          JfGrad,JndGrd,Work(ipMem1), nSO,Work(ipMem2),Mem2,
     &          Work(ipAux),nAux,Shijij)
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU2,TwoelWall2)
           Twoel2_CPU = Twoel2_CPU + TwoelCPU2-TwoelCPU1
           Twoel2_Wall = Twoel2_Wall + TwoelWall2-TwoelWall1
#endif
            If (iPrint.ge.15)
     &         Call PrGrad(' In Drvg1_2Center_RI: Grad',
     &                  Temp,nGrad,lIrrep,ChDisp,iPrint)
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
      Call GetMem('MemMax','Free','Real',ipMem1,MemMax)
      Call Free_Tsk(id)
      Call GetMem('ip_ij','Free','Inte',ip_ij,mij)
      Call GetMem('TMax','Free','Real',ipTMax,nTMax)
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
         Call GetMem('CijK','Free','Real',ip_CijK,2*lCijK)
         Call GetMem('A','Free','Real',ip_A,lA)
         If (iMP2Prpt.eq.2) Then
            Call GetMem('A_MP2(2)','Free','Real',ip_A_MP2(2),lA_MP2)
            Call GetMem('A_MP2(1)','Free','Real',ip_A_MP2(1),lA_MP2)
         End If
      End If
*
      Call Free_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Drvg1_2Center_RI')
      Return
      End
