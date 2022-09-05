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
* Copyright (C) 1990,1991,1993,1998,2005, Roland Lindh                 *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drv2El_2Center_RI(ThrAO,A_Diag,nSO_Aux,MaxCntr,SO2C)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals.                          *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. '93          *
*             Modified driver. Jan. '98                                *
*             Modified to 2-center ERIs for RI June '05                *
************************************************************************
      use Basis_Info, only: nBas_Aux
      use iSD_data
      use Wrj12, only: iOffA, Lu_A, SO2Ind
      use Index_arrays, only: iSO2Sh, nShBF
      use Gateway_Info, only: CutInt
      use RICD_Info, only: LDF
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External Integral_RI_2
#include "Molcas.fh"
#include "setup.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "nsd.fh"
#define _no_nShs_
#include "iTOffs.fh"
      Integer iAddr_AQ(0:7), kCol_Irrep(0:7)
      Logical Verbose, Indexation, FreeK2, DoGrad, DoFock
      Character Name_Q*6
      Real*8, Allocatable :: A_Diag(:)
      Integer, Allocatable:: SO2C(:)

      Real*8, Allocatable :: Tmp(:,:), TMax(:), TInt(:)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      Call StatusLine(' Seward:',' Computing 2-center RI integrals')
*                                                                      *
************************************************************************
*                                                                      *
*     Handle only the auxiliary basis set
*
      Call Set_Basis_Mode('Auxiliary')
      Call SetUp_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation.
*
      DoGrad=.False.
      DoFock=.False.
      Indexation = .True.
      Call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
*
      Call mma_Allocate(SO2Ind,nSOs,Label='SO2Ind')
      Call Mk_iSO2Ind(iSO2Sh,SO2Ind,nSOs,nSkal)
*
       nSO_Aux=nSOs-1
      If (LDF) Then
         Call mma_allocate(SO2C,nSO_Aux,Label='SO2C')
         MaxCntr=0
         Do i = 1, nSO_Aux
            iSh = iSO2Sh(i)
            iCenter=iSD(10,iSh)
            MaxCntr=Max(MaxCntr,iCenter)
            SO2C(i)=iCenter
         End Do
      Else
         MaxCntr=0
      End If
*
      nBfn2 = 0
      nBfnTot=0
      Do iIrrep = 0, nIrrep-1
         iTOffs(iIrrep+1) = nBfn2
*
         lJ=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) lJ=lJ-1
         nBfn2 = nBfn2 + lJ**2
         nBfnTot=nBfnTot+lJ
      End Do
      nA_Diag=nBfnTot
      Call mma_allocate(A_Diag,nA_Diag,Label='A_Diag')
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call mma_allocate(TMax,nSkal,Label='TMax')
      Call mma_allocate(Tmp,nSkal,nSkal,Label='Tmp')
      Call Shell_MxSchwz(nSkal,Tmp)

c     Call RecPrt('Tmp',' ',Tmp,nSkal,nSkal)

      TMax(:)=Tmp(:,nSkal)
      Call mma_deallocate(Tmp)
      TMax_all=Zero
      Do iS = 1, nSkal
         TMax_all=Max(TMax_all,TMax(iS))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Preallocate some core for Seward!
*
      Call mma_maxDBLE(MemSew)
      MemLow=Min(MemSew/2,1024*128)
      MemSew=Max(MemSew/10,MemLow)
      Call xSetMem_Ints(MemSew)
*                                                                      *
************************************************************************
*                                                                      *
*-----Temporary buffer for computed integrals, compute the largest
*     required buffer size and set up iOffA.
*
      nTInt=0
      Do jS = 1, nSkal-1
         nTInt = Max( nTInt,
     &                  nMemAm(nShBF,nIrrep,nSkal-1,jS,iOffA,
     &                         .True.) )
      End Do
      Call mma_allocate(TInt,nTInt,Label='TInt')
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*
*     Open files for the A-vectors, set iAddr_AQ, kCol_iIrrep and
*     iOffA(3,iIrrep).
*
      nBfnTot=0
      Do iIrrep = 0, nIrrep-1
         iOffA(3,iIrrep)=nBfnTot
         mB=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) mB = mB - 1
         nBfnTot=nBfnTot+mB
*
         iSeed=63+iIrrep
         Lu_A(iIrrep)=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'AVEC',iIrrep
         If (mB.ne.0) Call DaName_MF_WA(Lu_A(iIrrep),Name_Q)
*
         iAddr_AQ(iIrrep)=0
         kCol_Irrep(iIrrep)=0
      End Do
*
      iS = nSkal
      kS = nSkal
*
      Do jS = 1, nSkal-1
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        Initialize the buffer
*
         nTInt_=nMemAm(nShBF,nIrrep,nSkal-1,jS,iOffA,.True.)
         TInt(1:nTInt_)=Zero
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        First compute the A matrix
*
         Do lS = 1, jS
*
            Aint=TMax(jS)*TMax(lS)
            If (AInt.lt.CutInt) Go To 14
            Call Eval_IJKL(iS,jS,kS,lS,TInt,nTInt_,Integral_RI_2)
 14         Continue
*
*           Use a time slot to save the number of tasks and shell
*           quadrupltes process by an individual node
            Call SavStat(1,One,'+')
            Call SavStat(2,One,'+')
*
         End Do ! lS
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        Write the A-vectors to disk
*
         Do iIrrep = 0, nIrrep-1
            mB = iOffA(2,iIrrep)           ! # of bf of shell jS
            If (mB.ne.0) Then
*
               ip_A_n = 1 + iOffA(1,iIrrep)
               iAddr=iAddr_AQ(iIrrep) ! Disk address
*
               nB = nBas_Aux(iIrrep)
               If (iIrrep.eq.0) nB = nB - 1 ! subtract dummy af
               Do kCol = 1+kCol_Irrep(iIrrep), mB+kCol_Irrep(iIrrep)
*
*                 Write the A-vector to file
*
                  Call dDaFile(Lu_A(iIrrep),1,TInt(ip_A_n),kCol,iAddr)

                  ipAs_Diag=1+iOffA(3,iIrrep)+kCol-1
                  A_Diag(ipAs_Diag)=TInt(ip_A_n+kCol-1)
                  nZero=nB-kCol
                  If (nZero.ne.0) Call dDaFile(Lu_A(iIrrep),0,
     &                                         TInt(ip_A_n),
     &                                         nZero,iAddr)
*
                  ip_A_n = ip_A_n + kCol
               End Do
*
               kCol_Irrep(iIrrep)=kCol_Irrep(iIrrep)+mB
               iAddr_AQ(iIrrep)=iAddr
            End If
*
         End Do ! iIrrep
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
      End Do    ! jS
*----------------------------------------------------------------------*
*                                                                      *
*     Release the Seward core memory, the buffer, etc.
*
      Call Free_iSD()
      Call xRlsMem_Ints
      Call mma_deallocate(TInt)
      Call mma_deallocate(TMax)
      Call mma_deallocate(SO2Ind)
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate integral environment.
*
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
      Call SavTim(1,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
