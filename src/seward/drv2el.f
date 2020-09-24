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
* Copyright (C) 1990,1991,1993,1998, Roland Lindh                      *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drv2El(Integral_WrOut,ThrAO)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals.                          *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              Timing                                                  *
*              Setup_Ints                                              *
*              Eval_Ints                                               *
*              Term_Ints                                               *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Modified for k2 loop. August '91                         *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. '93          *
*             Modified driver. Jan. '98                                *
************************************************************************
      use iSD_data
      use Basis_Info, only: dbsc
      use Real_Info, only: CutInt
      Implicit Real*8 (A-H,O-Z)
      External Integral_WrOut, Rsv_GTList
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
      Parameter(nTInt=1,mDens=1)
      Real*8 Dens(mDens), Fock(mDens), TInt(nTInt)
      Integer iTOffs(8,8,8)
      Logical Verbose, Indexation, FreeK2,
     &        W2Disc, PreSch, DoIntegrals, DoFock, DoGrad,
     &        FckNoClmb, FckNoExch, Rsv_GTList, Triangular
      Character*72 SLine
      Real*8, Dimension(:,:), Allocatable :: TMax
      Integer, Dimension(:,:), Allocatable :: Pair_Index
      Dimension ExFac(1),FckNoClmb(1),FckNoExch(1)
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      Call QEnter('Drv2El')
      SLine='Computing 2-electron integrals'
      Call StatusLine(' Seward:',SLine)
*                                                                      *
************************************************************************
*                                                                      *
      ExFac=One
      Nr_Dens=1
      DoIntegrals=.True.
      DoFock=.False.
      DoGrad=.False.
      FckNoClmb=.False.
      FckNoExch=.False.
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do not generate
*     tables for indexation.
*
      Indexation = .False.
      Call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
*                                                                      *
************************************************************************
*                                                                      *
      Thize=Zero               ! Not used for conventional integrals
      PreSch=.True.            ! Not used for conventional integrals
      Disc_Mx=Zero             ! Not used for conventional integrals
*
      Disc=Zero
      Dix_Mx=Zero
      TskHi=Zero
      TskLw=Zero
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call mma_allocate(TMax,nSkal,nSkal)
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
      Call mma_allocate(Pair_Index,2,nSkal*(nSkal+1)/2)
      nij=0
      Do iS = 1, nSkal
         Do jS = 1, iS
            If (TMax_All*TMax(iS,jS).ge.CutInt) Then
               nij = nij + 1
               Pair_Index(1,nij)=iS
               Pair_Index(2,nij)=jS
            End If
         End Do
      End Do
      P_Eff=DBLE(nij)
*                                                                      *
************************************************************************
*                                                                      *
      Triangular=.True.
      Call Alloc_TList(Triangular,P_Eff)
      Call Init_TList(Triangular,P_Eff)
      Call Init_PPList
      Call Init_GTList
      iOpt=0
*
      PP_Eff=P_Eff**2
      PP_Eff_delta=0.10D0*PP_Eff
      PP_Count=Zero
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*
*     big loop over individual tasks distributed over individual nodes
*
 10   Continue
*     make reservations of a tesk in global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_GTlist(TskLw,TskHi,iOpt,W2Disc)) Go To 11
      W2Disc=.False.
*
*     Now do a quadruple loop over shells
*
      ijS = Int((One+sqrt(Eight*TskLw-Three))/Two)
      iS = Pair_Index(1,ijS)
      jS = Pair_Index(2,ijS)
      klS = Int(TskLw-DBLE(ijS)*(DBLE(ijS)-One)/Two)
      kS = Pair_Index(1,klS)
      lS = Pair_Index(2,klS)
      Count=TskLw
*
      If (Count-TskHi.gt.1.0D-10) Go To 12
 13   Continue
*
*     Logic to avoid computing integrals in a mixed muonic and
*     electronic basis.
*
      iCnttp=iSD(13,iS)
      jCnttp=iSD(13,jS)
      If (dbsc(iCnttp)%fMass.ne.dbsc(jCnttp)%fMass) Go To 14
      kCnttp=iSD(13,kS)
      lCnttp=iSD(13,lS)
      If (dbsc(kCnttp)%fMass.ne.dbsc(lCnttp)%fMass) Go To 14
*
      S_Eff=DBLE(ijS)
      T_Eff=DBLE(klS)
      ST_Eff=S_Eff*(S_Eff-One)/2D0 + T_Eff
      If (ST_Eff.ge.PP_Count) Then
         Write (SLine,'(A,F5.2,A)') 'Computing 2-electron integrals,',
     &        ST_Eff/PP_Eff*100d0,'% done so far.'
         Call StatusLine(' Seward:',SLine)
         PP_Count = PP_Count + PP_Eff_delta
      End If
*
*
         Aint=TMax(iS,jS)*TMax(kS,lS)
         If (AInt.lt.CutInt) Go To 14
         Call Eval_Ints_New_Internal
     &                  (iS,jS,kS,lS,TInt,nTInt,
     &                   iTOffs,Integral_WrOut,
* the following are dummy arguments
     &                   Dens,Fock,mDens,ExFac,Nr_Dens,
     &                   FckNoClmb,FckNoExch,
     &                   Thize,W2Disc,PreSch,Dix_Mx,Disc,
     &                   Count,DoIntegrals,DoFock)
 14      Continue
         Count=Count+One
         If (Count-TskHi.gt.1.0D-10) Go To 12
         klS = klS + 1
         If (klS.gt.ijS) Then
            ijS = ijS + 1
            klS = 1
         End If
         iS = Pair_Index(1,ijS)
         jS = Pair_Index(2,ijS)
         kS = Pair_Index(1,klS)
         lS = Pair_Index(2,klS)
         Go To 13
*
*     Task endpoint
 12   Continue
*
*     Use a time slot to save the number of tasks and shell
*     quadrupltes process by an individual node
      Call SavStat(1,One,'+')
      Call SavStat(2,TskHi-TskLw+One,'+')
      Go To 10
 11   Continue
*     End of big task loop
      Call CWTime(TCpu2,TWall2)
      Call SavTim(1,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_GTList
      Call Free_PPList
      Call Free_TList
*
      Call mma_deallocate(Pair_Index)
      Call mma_deallocate(TMax)
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate integral environment.
*
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
      Call Free_iSD()
      Call QExit('Drv2El')
      Return
      End
