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
      SubRoutine Drv2El_2Center_RI(Integral_WrOut,ThrAO,ipA_Diag,
     &                             nSO_Aux,MaxCntr,ipSO2C)
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
      use Wrj12
      use Index_arrays, only: iSO2Sh, nShBF
      use Real_Info, only: CutInt
      use RICD_Info, only: LDF
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External Integral_WrOut
#include "Molcas.fh"
#include "setup.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "nsd.fh"
#define _no_nShs_
#include "iTOffs.fh"
      Integer iAddr_AQ(0:7), kCol_Irrep(0:7)
      Logical Verbose, Indexation, FreeK2, DoGrad, DoFock
      Character Name_Q*6
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*----- Statement functions
*
      TMax(i)=Work(ipTMax-1+i)
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
         Call GetMem('SO2C','Allo','Inte',ipSO2C,nSO_Aux)
         MaxCntr=0
         Do i = 1, nSO_Aux
            iSh = iSO2Sh(i)
            iCenter=iSD(10,iSh)
            MaxCntr=Max(MaxCntr,iCenter)
            iWork(ipSO2C+i-1)=iCenter
         End Do
      Else
         MaxCntr=0
         ipSO2C=ip_Dummy
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
      Call GetMem('A_Diag','Allo','Real',ipA_Diag,nA_Diag)
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call GetMem('TMax','Allo','Real',ipTMax,nSkal)
      Call Allocate_Work(ip_Tmp,nSkal**2)
      Call Shell_MxSchwz(nSkal,Work(ip_Tmp))

c      Call RecPrt('ip_Tmp',' ',Work(ip_Tmp),nSkal,nSkal)

      call dcopy_(nSkal,Work(ip_Tmp+(nSkal-1)*nSkal),1,Work(ipTMax),1)
      Call Free_Work(ip_Tmp)
      TMax_all=Zero
      Do iS = 1, nSkal
         TMax_all=Max(TMax_all,TMax(iS))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Preallocate some core for Seward!
*
      Call GetMem('MaxMem','Max','Real',iDummy,MemSew)
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
      Call GetMem('Am','Allo','Real',ipTInt,nTInt)
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
         Call FZero(Work(ipTInt),nTInt_)
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        First compute the A matrix
*
         Do lS = 1, jS
*
            Aint=TMax(jS)*TMax(lS)
            If (AInt.lt.CutInt) Go To 14
            Call Eval_IJKL(iS,jS,kS,lS,Work(ipTInt),nTInt_,
     &                    Integral_WrOut)
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
               ip_A_n=ipTInt + iOffA(1,iIrrep)
               iAddr=iAddr_AQ(iIrrep) ! Disk address
*
               nB = nBas_Aux(iIrrep)
               If (iIrrep.eq.0) nB = nB - 1 ! subtract dummy af
               Do kCol = 1+kCol_Irrep(iIrrep), mB+kCol_Irrep(iIrrep)
*
*                 Write the A-vector to file
*
                  Call dDaFile(Lu_A(iIrrep),1,Work(ip_A_n),
     &                         kCol,iAddr)

                  ipAs_Diag=ipA_Diag+iOffA(3,iIrrep)+kCol-1
                  Work(ipAs_Diag)=Work(ip_A_n+kCol-1)
                  ipAs_Diag=ipAs_Diag+1
                  nZero=nB-kCol
                  If (nZero.ne.0) Call dDaFile(Lu_A(iIrrep),0,
     &                                         Work(ip_A_n),
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
      Call GetMem('Am', 'Free','Real',ipTInt,nTInt)
      Call GetMem('TMax','Free','Real',ipTMax,nSkal)
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
