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
* Copyright (C) Ben Swerts                                             *
************************************************************************
CStart Molcas
      SubRoutine Drvg_FAIEMP(Grad,Temp,nGrad)
celse
c;      SubRoutine Drvg_FAIEMP(Grad,Temp,nGrad,fock,fock1,d1ao,d1ao1)
cend
************************************************************************
*                                                                      *
*  Object: driver for the derivatives of central-fragment              *
*          two-electron integrals. The four outermost loops            *
*          will controll the type of the two-electron integral, eg.    *
*          (ss|ss), (sd|pp), etc. The next four loops will generate    *
*          list of symmetry distinct centers that do have basis func-  *
*          tions of the requested type.                                *
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
*     Author: Ben Swerts                                               *
*                                                                      *
*     based on Drvg1                                                   *
*                                                                      *
************************************************************************
      use k2_setup
      use iSD_data
      use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr

      Implicit None
      External King, Rsv_GTList, MPP
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
*     Local arrays
      Integer  nGrad
      Real*8   Coor(3,4), Grad(nGrad), Temp(nGrad)
      Integer  iAnga(4), iCmpa(4), iShela(4),iShlla(4),
     &         iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4)
      Integer  nHrrTb(0:iTabMx,0:iTabMx,2)
      Logical  EQ, Shijij, AeqB, CeqD, lDummy,
     &         DoGrad, DoFock, Indexation,
     &         JfGrad(3,4), ABCDeq, No_Batch, King, Rsv_GTList, MPP,
     &         FreeK2, Verbose, Triangular
      Character*7 Format*72
      Logical  lNoSkip
      Integer  nBas_Valence(0:7)
*
      Integer  iSD4(0:nSD,4)
      Integer  MemMax,MemPrm
      save     MemPrm
*
      Integer  iRout,iPrint,nBT,nBVT,idum,idum1,i,j,iAng,iBasi,iBasn
      Integer  iS,jS,iBasAO,iBsInc,iCar,ijklA,ijS,Indij,iOpt,ijMax
      Integer  ip_ij,ipEI,ipEta,ipiEta,ipCffi,ipMem1,ipMem2,ipP,ipQ
      Integer  iPrem,iPren,ipxA,ipxB,ipxG,ipxD,ipZi,Mem1,Mem2,iPrimi
      Integer  iPrInc,ipTMax,jAng,iSh,jBasAO,jBasj,jBasn,jBsInc,jpCffj
      Integer  jPrInc,k2ij,k2kl,jPrimj,kBasAO,kBasn,kBask,kBsInc
      Integer  kBtch,klS,kpCffk,kPrimk,kPrInc,kS,lBasAO,lBasl,lBasn
      Integer  lBsInc,lpCffl,lPriml,lPrInc,mBtch,lS,mdci,mdcj,mdck,mdcl
      Integer  MemPSO,nab,ncd,nDCRR,nDCRS,nEta,nHmab,nHmcd,nHrrab
      Integer  nij,nijkl,nPairs,nQuad,nRys,nSkal,nSkal_Fragments
      Integer  nSkal_Valence,nSO,nZeta,nBtch
      Real*8   TMax,PMax,ExFac,CoulFac,Aint,Count,P_Eff,Prem,Pren
      Real*8   TCpu1,TCpu2,ThrAO,TMax_all,TskHi,TskLw,TWall1,TWall2
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      TMax(i,j)=Work((j-1)*nSkal+i+ipTMax-1)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 203
      iPrint = nPrint(iRout)
      iFnc(1)=0
      iFnc(2)=0
      iFnc(3)=0
      iFnc(4)=0
      PMax=Zero
      idum=0
      idum1=0
      ExFac=One
      CoulFac=One
      Call QEnter('Drvg_FAIEMP')
*
*     Handle both the valence and the fragment basis set
*
      Call Set_Basis_Mode('Valence')
      Call Nr_Shells(nSkal_Valence)
      Call Set_Basis_Mode('WithFragments')
      Call SetUp_iSD
      nBT = 0
      nBVT = 0
      do i = 0, nIrrep - 1
        nBas_Valence(i) = nBas(i)
        nBVT = nBVT + nBas(i)*(nBas(i)+1)/2
        nBas(i) = nBas(i) + nBas_Frag(i)
        nBT = nBT + nBas(i)*(nBas(i)+1)/2
      enddo
*                                                                      *
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
      nSkal_Fragments=nSkal-nSkal_Valence
      If(iPrint.ge.99) write(6,*) 'nSkal, nSkal_Valence, nSkal_Frag = ',
     &                             nSkal, nSkal_Valence, nSkal_Fragments
*                                                                      *
************************************************************************
*                                                                      *
*-----Prepare handling of the combined two-particle density.
*
      Call PrepP_FAIEMP(nBas_Valence, nBT, nBVT)
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
      Call GetMem('TMax','Allo','Real',ipTMax,nSkal**2)
      Call Shell_MxSchwz(nSkal,Work(ipTMax))
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
      Call GetMem('ip_ij','Allo','Inte',ip_ij,nSkal*(nSkal+1))
      nij=0
      Do iS = 1, nSkal
         Do jS = 1, iS
            If (TMax_All*TMax(iS,jS).ge.CutInt) Then
               nij = nij + 1
               iWork((nij-1)*2+ip_ij  )=iS
               iWork((nij-1)*2+ip_ij+1)=jS
            End If
         End Do
      End Do
      P_Eff=DBLE(nij)
*                                                                      *
************************************************************************
*                                                                      *
*-------Compute FLOP's for the transfer equation.
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
      call dcopy_(nGrad,[Zero],0,Temp,1)
      If (iPrint.ge.15) Call PrGrad(' In Drvg_FAIEMP: Total Grad (1)',
     &                              Grad,nGrad,lIrrep,ChDisp,iprint)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_MaxDBLE(MemMax)
      Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
      ipMem1=1
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
      iS = iWork((ijS-1)*2+ip_ij)
      jS = iWork((ijS-1)*2+ip_ij+1)
      klS = Int(TskLw-DBLE(ijS)*(DBLE(ijS)-One)/Two)
      kS = iWork((klS-1)*2+ip_ij)
      lS = iWork((klS-1)*2+ip_ij+1)
      Count=TskLw
      Call CWTime(TCpu1,TWall1)
  13  Continue
         If (Count-TskHi.gt.1.0D-10) Go To 12
*
         Aint=TMax(iS,jS)*TMax(kS,lS)

         lNoSkip = AInt.ge.CutInt
* only calculate needed integrals and only update the valence part of the
* Fock matrix (iS > nSkal_Valence, lS <= nSkal_Valence, jS and kS
* belonging to different regions)
         If(jS.le.nSkal_Valence) Then
           lNoSkip = lNoSkip.and.kS.gt.nSkal_Valence
         Else
           lNoSkip = lNoSkip.and.kS.le.nSkal_Valence
         End If
         lNoSkip = lNoSkip.and.lS.le.nSkal_Valence
         if(.Not.lNoSkip) Go To 14
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
         if(nIrrep.eq.1.and.ABCDeq.and.Mod(ijklA,2).eq.1
     &     .and.iPrint.gt.15) write(6,*) 'ABCDeq & ijklA'
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
           Call PGet0(iCmpa,iShela,
     &                iBasn,jBasn,kBasn,lBasn,Shijij,
     &                iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,
     &                iFnc(1)*iBasn,iFnc(2)*jBasn,
     &                iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,
     &                Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
           If (AInt*PMax.lt.CutInt) Go To 430
*
*----------Compute gradients of shell quadruplet
*
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
     &          Mem_DBLE(ipZeta),Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,
     &          Mem_DBLE(ipEta), Mem_DBLE(ipEI),Mem_DBLE(ipQ),nEta,
     &          Mem_DBLE(ipxA),Mem_DBLE(ipxB),
     &          Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,
     &          JfGrad,JndGrd,Sew_Scr(ipMem1), nSO,Sew_Scr(ipMem2),Mem2,
     &          Aux,nAux,Shijij)
*
*
            If (iPrint.ge.15)
     &         Call PrGrad(' In Drvg_FAIEMP: Grad',
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
 14      Continue
         Count=Count+One
         klS = klS + 1
         If (klS.gt.ijS) Then
            ijS = ijS + 1
            klS = 1
         End If
         iS = iWork((ijS-1)*2+ip_ij  )
         jS = iWork((ijS-1)*2+ip_ij+1)
         kS = iWork((klS-1)*2+ip_ij  )
         lS = iWork((klS-1)*2+ip_ij+1)
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
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Sew_Scr)
      Call Free_GTList
      Call Free_PPList
      Call Free_TList
      Call GetMem('ip_ij','Free','Inte',ip_ij,nSkal*(nSkal+1))
      Call GetMem('TMax','Free','Real',ipTMax,nSkal**2)
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
     &   ' A total of', Pren,' entities were prescreened and',
     &                  Prem,' were kept.'
      End If
* Accumulate the final results
      Call DScal_(nGrad,Half,Temp,1)
      If(iPrint.ge.15) Call PrGrad('The FAIEMP 2-electron Contribution',
     &                             Temp,nGrad,lIrrep,ChDisp,iPrint)
      call daxpy_(nGrad,One,Temp,1,Grad,1)
*
      Call Free_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Drvg_FAIEMP')
      Return
      End
