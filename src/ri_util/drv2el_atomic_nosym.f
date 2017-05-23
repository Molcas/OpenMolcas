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
      SubRoutine Drv2El_Atomic_NoSym(Integral_WrOut,ThrAO,iCnttp,jCnttp,
     &                               TInt,nTInt,
     &                               In_Core,ADiag,LuA,ijS_req,
     &                               Keep_Shell)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals.                          *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              Timing                                                  *
*              Setup_Ints                                              *
*              Eval_IJKL                                               *
*              Term_Ints                                               *
*              QExit                                                   *
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
*                                                                      *
************************************************************************
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
      External Integral_WrOut
#include "itmax.fh"
#include "info.fh"
#include "lundio.fh"
#include "nsd.fh"
#include "setup.fh"
#include "print.fh"
#include "real.fh"
#include "shinf.fh"
#include "wrj12.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
#define _no_nShs_
#include "iTOffs.fh"
      Real*8, Allocatable :: TInt(:), ADiag(:)
      Logical Verbose, Indexation, FreeK2, DoGrad, DoFock,
     &        In_Core, Out_of_Core, Only_DB, Do_RI_Basis, Do_ERIs
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      Call QEnter('Drv2El_Atomic_NoSym')
*                                                                      *
************************************************************************
*                                                                      *
*     Temporary modifications to facilitate atomic calculations
*
      nIrrep_Save=nIrrep
      nIrrep=1
      Petite=.True.
      iWROpt_Save=iWROpt
      iWROpt=1
*
      Do_RI_Basis = AuxCnttp(iCnttp)
*
      Call Set_Basis_Mode_Atomic(iCnttp,jCnttp)
      Call Setup_iSD()
*
      If (Do_RI_Basis .and. ijS_req.eq.0) Then
         Call WarningMessage(2,
     &               'Do_RI_Basis .and. ijS_req.eq.0')
         Call Abend()
      End If
C     Write (6,*) 'Do_RI_Basis=',Do_RI_Basis
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do not generate
*     tables for indexation.
*
      DoGrad=.False.
      DoFock=.False.
      Indexation = .False.
      Call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
C     Write (6,*) 'nSkal=',nSkal
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of pairs
*
      Call GetMem('ip_ij','Allo','Inte',ip_ij,nSkal*(nSkal+1))
      nij=0
      nBfn=0
      If (Do_RI_Basis) Then
         iS=nSkal   ! Dummy shell
         Do jS = 1, nSkal-1
            nij = nij +1
            iWork((nij-1)*2+ip_ij  )=iS
            iWork((nij-1)*2+ip_ij+1)=jS
            nBfn=nBfn+iSD(2,jS)*iSD(3,jS)
         End Do
      Else
         Do iS = 1, nSkal
            Do jS = 1, iS
               nij = nij +1
               iWork((nij-1)*2+ip_ij  )=iS
               iWork((nij-1)*2+ip_ij+1)=jS
            End Do
         End Do
      End If
C     Write (6,*) 'nij=',nij
*                                                                      *
************************************************************************
*                                                                      *
*     Preallocate some core for Seward!
*
      Call GetMem('MaxMem','Max','Real',iDummy,MemSew)
*
      MemLow=Min(MemSew/2,1024*128)
      MemSew=Max(MemSew/10,MemLow)
      Call xSetMem_Ints(MemSew)
*                                                                      *
************************************************************************
*                                                                      *
*     Determine if only diagonal block should be computed.
*     This option is forced when the auxiliary basis set is transformed
*     to the Cholesky basis!
*
      Only_DB=ijS_req.ne.0 .or. Do_RI_Basis
*                                                                      *
************************************************************************
*                                                                      *
*     Choose between in-core and out-of-core options
*
      Call GetMem('MaxMem','Max','Real',iDummy,MemT)
      MemT=MemT/2
*
      If (Only_DB) Then
*
*        Only diagonal block
*
         nTInt=0
         If (Do_RI_Basis) Then
            Do iS = 1, nSkal-1 ! Skip the dummy shell
               nBfn_i=iSD(2,iS)*iSD(3,iS)
               If (iS.eq.ijS_req) nTInt=nBfn_i
            End Do
            If (nTInt.eq.0) Then
               Call WarningMessage(2,
     &                     'Drv2el_atomic_nosym: nTInt.eq.0')
               Call Abend()
            End If
*
            Call GetMem('SO2Ind','Allo','Inte',ipSO2Ind,nBfn)
            Do iBfn = 1, nBfn
               iWork(ipSO2Ind+iBfn-1)=iBfn
            End Do
            nSOs=nBfn
         Else
            Do iS = 1, nSkal
               nBfn_i=iSD(2,iS)*iSD(3,iS)
               Do jS = 1, iS-1
                  ijS = iS*(iS-1)/2 + jS
                  nBfn_j=iSD(2,jS)*iSD(3,jS)
                  If (ijS.eq.ijS_req) nTInt=nBfn_i*nBfn_j
               End Do
               ijS = iS*(iS+1)/2
               If (ijS.eq.ijS_req) nTInt=nBfn_i*(nBfn_i+1)/2
            End Do
         End If
         mTInt=nTInt
         nTInt2=mTInt*nTInt
         If (nTInt2.gt.MemT) Then
            Call WarningMessage(2,'Not enough memory!')
            Call Abend()
         End If
*
         iTOffs(1)=0
*
         In_Core=.True.       ! no out-of-core option needed.
         Out_of_core=.False.
         mTInt2=nTInt2
*
      Else
*
*        All blocks
*
         nTInt=nBas(0)*(nBas(0)+1)/2
         mTInt=nTInt
         nTInt2=nTInt**2
         In_Core=nTInt2.le.MemT
         If (Force_out_of_Core) In_Core=.False.
         Out_of_Core=.NOT.In_Core
*
*        Compute the size of the array, TInt, to write the integrals to.
*
         If (Out_of_Core) Then
*
*           Find the larges block for a fixed shell pair.
*
            mTInt=0
            Do iS = 1, nSkal
               nBfn_i=iSD(2,iS)*iSD(3,iS)
               Do jS = 1, iS-1
                  nBfn_j=iSD(2,jS)*iSD(3,jS)
                  mTInt=Max(mTInt,nBfn_i*nBfn_j)
               End Do
               mTInt=Max(mTInt,nBfn_i*(nBfn_i+1)/2)
            End Do
            nTInt2=mTInt*nTInt
*
*           Open file for the A-vectors
*
            iAddr=0
            iSeed=63
            LuA=IsFreeUnit(iSeed)
            Call DaName_MF_WA(LuA,'AVEC0')
*
            Call mma_allocate(ADiag,nTInt,label='ADiag')
            Call FZero(ADiag,nTInt)
*
         Else
*
*           In-core option
*
            mTInt2=nTInt2
            iTOffs(1)=0  ! Offset permanently set to zero!
*
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      iTOffs(2)=nTInt  ! # of rows in TInt
      iTOffs(3)=mTInt  ! # of colums in TInt
      Call mma_allocate(TInt,nTInt2,label='TInt')
      If (In_Core) Call FZero(TInt,nTInt2)
*                                                                      *
************************************************************************
*                                                                      *
*     Now do a quadruple loop over shells
*
      iTOff=0
      iTOffs(4)=0    ! Offset to the ij set
      Do ijS = 1, nij
         iS = iWork((ijS-1)*2+ip_ij)
         jS = iWork((ijS-1)*2+ip_ij+1)
*
         nBfn_i=iSD(2,iS)*iSD(3,iS)
         nBfn_j=iSD(2,jS)*iSD(3,jS)
         ijAng =iSD(1,iS)+iSD(1,jS)
*
         If (Out_of_Core) Then
            If (iS.eq.jS) Then
               mTInt=nBfn_i*(nBfn_i+1)/2
            Else
               mTInt=nBfn_i*nBfn_j
            End If
            mTInt2=mTInt*nTInt
            Call FZero(TInt,mTInt2)
            iTOffs(1)=iTOff
            iTOffs(3)=mTInt
         End If
*
         iTOffs(5)=0    ! Offset to the kl set
         Do klS = 1, ijS
            kS = iWork((klS-1)*2+ip_ij)
            lS = iWork((klS-1)*2+ip_ij+1)
*
            nBfn_k=iSD(2,kS)*iSD(3,kS)
            nBfn_l=iSD(2,lS)*iSD(3,lS)
            klAng =iSD(1,kS)+iSD(1,lS)
*
*           For Only_DB compute the shell quadruplet (ijS_req|ijS_req)
*
            Do_ERIs = .NOT.Only_DB .or.
     &                (ijS.eq.ijS_req.and.klS.eq.ijS_req)
*
*           If high angular combination of the product basis is skipped
*           do not compute the contributions.
*
            Do_ERIs = Do_ERIs .and. ijAng.le.Keep_Shell .and.
     &                              klAng.le.Keep_Shell
*
            If (Do_ERIs) Then
               Call Eval_IJKL(iS,jS,kS,lS,TInt,mTInt2,
     &                       Integral_WrOut)
            End If
*
            If (.NOT.Only_DB) Then
               If (kS.eq.lS) Then
                  iTOffs(5)=iTOffs(5)+nBfn_k*(nBfn_k+1)/2
               Else
                  iTOffs(5)=iTOffs(5)+nBfn_k*nBfn_l
               End If
            End If
*
         End Do
*
*        For out-of-core version write the integrals to disk!
*        Pick up the diagonal elements.
*
         If (Out_of_Core) Then
            Call dDaFile(LuA,1,TInt,mTInt2,iAddr)
            call dcopy_(mTInt,TInt(1+iTOff),nTInt+1,
     &                  ADiag(1+iTOff),1)
            iTOff=iTOff+mTInt
         End If
*
         If (.NOT.Only_DB) Then
            If (iS.eq.jS) Then
               iTOffs(4)=iTOffs(4)+nBfn_i*(nBfn_i+1)/2
            Else
               iTOffs(4)=iTOffs(4)+nBfn_i*nBfn_j
            End If
         End If
*
      End Do      !    ijS
*
      If (Do_RI_Basis) Call Free_iWork(ipSO2Ind)
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*
      If (Out_of_Core) Call mma_deallocate(TInt)
      Call xRlsMem_Ints
      Call GetMem('ip_ij','Free','Inte',ip_ij,nSkal*(nSkal+1))
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
*     Square TInt from upper triangular to full.
*
      If (In_Core.and..Not.Do_RI_Basis) Then
*
         Do iTInt = 1, nTInt
            Do jTInt = 1, iTInt-1
               ij = (jTInt-1)*nTInt + iTInt
               ji = (iTInt-1)*nTInt + jTInt
               TInt(ij)=TInt(ji)
            End Do
         End Do
C        Call RecPrt('Drv2El_atomic: TInt',' ',
C    &               TInt,nTInt,nTInt)
*
      Else If (.Not.Do_RI_Basis) Then
*
         nij=nBas(0)*(nBas(0)+1)/2
         Call GetMem('MemMax','Max','Real',iDummy,MaxMem)
         Call Square_A(LuA,nij,MaxMem,Force_Out_of_Core)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      nIrrep=nIrrep_Save
      Petite = nIrrep.eq.1
      iWROpt=iWROpt_Save
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Drv2El_Atomic_NoSym')
      Return
      End
