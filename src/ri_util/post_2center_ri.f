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
      SubRoutine Post_2Center_RI(ipA_Diag)
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
      use Wrj12
      use Temporary_Parameters, only: force_out_of_core
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "setup.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer  nDmA(0:7),  nDmB(0:7)
      Logical Out_of_Core
      Character Name_Q*6

      Real*8, Allocatable :: Scr(:), X(:), Z(:)
      Integer, Allocatable :: iDiag(:)
      Real*8, Allocatable, Target :: Am(:), Qm(:), A_k(:), Q_k(:)
      Real*8, Pointer :: A_l(:)=>Null(), Q_l(:)=>Null()
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
      nScr=0
      nBfn2 = 0
      nBfnTot=0
      Do iIrrep = 0, nIrrep-1
         lJ=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) lJ=lJ-1
         nDmA(iIrrep)=lJ
         nDmB(iIrrep)=0
         nScr=Max(nScr,3*lJ)
         nBfn2 = nBfn2 + lJ**2
         nBfnTot=nBfnTot+lJ
      End Do
      nA_Diag=nBfnTot
*
      Call mma_maxDBLE(MaxMem)
*                                                                      *
************************************************************************
*                                                                      *
*     Fill in the lower part of the A matrix as it is stored on disk.
*
      Do iIrrep = 0, nIrrep-1
         nB = nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB = nB - 1 ! subtract dummy af
         Call Square_A(Lu_A(iIrrep),nB,MaxMem,Force_Out_of_Core)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Pivoting of the A matrix
*
      Call mma_allocate(iDiag,nA_Diag,Label='iDiag')
      Call mma_maxDBLE(MaxMem2)
*
      If (Force_Out_of_Core) MaxMem2=3*nBfnTot
*     lScr=Min(MaxMem2,nScr)
      lScr=Max(MaxMem2-(nScr/3),nScr)
      Call mma_allocate(Scr,lScr,Label='Scr')
*
      Call SORT_mat(irc,Work(ipA_Diag),nDmA,nDmB,iDiag,nIrrep,
     &                  LU_A,'GePivot',lScr,Scr)
      ichk=0
      Do iIrrep = 0, nIrrep-1
         nChV(iIrrep)=nDmB(iIrrep)
         ichk=ichk+Min(1,nDmA(iIrrep)-nDmB(iIrrep))
      End Do
      If (ichk.ne.0) Then
         write(6,*)
         write(6,*)'Post_2Center_RI'
         write(6,*)'Detected lin. dependences in the auxiliary basis.'
         Write(6,'(A,8I6)')
     & ' # of AuxBas before l. d. removal: ',(nDmA(i),i=0,nIrrep-1)
         Write(6,'(A,8I6)')
     & ' # of AuxBas after  l. d. removal: ',(nDmB(i),i=0,nIrrep-1)
         write(6,*)
      EndIf
*
      Call SORT_mat(irc,Work(ipA_Diag),nDmA,nDmB,iDiag,nIrrep,
     &                  LU_A,'DoPivot',lScr,Scr)
*
*     Note: after the 'DoPivot' call to Sort_mat, the A-matrix is
*           no longer stored as squared but as upper-triangular
*
      Call mma_deallocate(Scr)
      Call GetMem('A_Diag','Free','Real',ipA_Diag,nA_Diag)
      ipA_Diag=ip_Dummy  ! Dummy memory pointer.
*
************************************************************************
*     A-vectors are now on disk. Go ahead and compute the Q-vectors!
************************************************************************
*
      ThrQ=1.0d-14 ! Threshold for Inv_Cho_Factor
*
      Do iIrrep = 0, nIrrep-1
c         nB=nBas_Aux(iIrrep)
c         If (iIrrep.eq.0) nB = nB - 1
         nB=nDmB(iIrrep)
         If (nB.eq.0) Go To 777
         nQm=nB*(nB+1)/2
*
         nXZ=nB
         nQm_full= nB*(nB+1)/2
*
         If (Force_Out_of_Core) MaxMem=(8*(2*nQm_full+5*nXZ))/10
         Out_of_Core=2*nQm_full+5*nXZ.gt.MaxMem
*
         If (Out_Of_Core) Then
            mQm=(nQm*MaxMem-5*nXZ)/(2*nQm_full)
            a=One
            b=-Two*DBLE(mQm)
            mB=INT(-a/Two + Sqrt( (a/Two)**2 - b ))
            kQm=mB*(mB+1)/2
            If (kQm.gt.mQm) Then
               Call WarningMessage(2,'Error in Post_2Center_RI')
               Write (6,*) 'kQm.gt.mQm!'
               Write (6,*) 'MaxMem=',MaxMem
               Write (6,*) 'nQm,mQm,kQm=',nQm,mQm,kQm
               Write (6,*) 'nB,mB=',nB,mB
               Call Abend()
            End If
         Else
            mB = nB
            kQm = nQm
         End If
*
         lQm=kQm
         lAm=lQm
*
         If (lQm.lt.1) Then
            Call WarningMessage(2,'Error in Post_2Center_RI')
            Write (6,*) 'lQm.lt.1'
            Call Abend()
         End If
*
*        Some of memory for scratch arrays for Inv_Cho_Factor
*        Allocate memory for the A- and Q-vectors and initialize.
*
         lScr=nXZ
         Call mma_allocate(Scr,lScr,Label='Scr')
         Call mma_allocate(Z,nXZ,Label='Z')
         Call mma_allocate(X,nXZ,Label='X')
         Call mma_allocate(Am,lAm,Label='Am')
         ip_Am = ip_of_Work(Am(1))
         Call mma_allocate(Qm,lQm,Label='Qm')
         ip_Qm = ip_of_Work(Qm(1))
         Call mma_allocate(A_k,nXZ,Label='A_k')
         Call mma_allocate(Q_k,nXZ,Label='Q_k')
         ip_Q_k = ip_of_Work(Q_k(1))
         ip_A_k = ip_of_Work(A_k(1))
*
         Am(:)=Zero
         Qm(:)=Zero
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        Process the A_ks to generate Q_ks.
*
         iSeed=55+iIrrep
         Lu_Q(iIrrep)=IsFreeUnit(iSeed)
         Write(Name_Q,'(A4,I2.2)') 'QMAT',iIrrep
         Call DaName_MF_WA(Lu_Q(iIrrep),Name_Q)
*
         iAddr=0
         nMem=mB
         Do kCol = 1, nB
*
            iAddr_=iAddr
            If (kCol.le.nMem) Then
*              Point to A_k in Am
               iOff = (kCol-1)*kCol/2
               ip_A_l=ip_Am+iOff
               If (kCol.eq.1) Then
                  nAm=nMem*(nMem+1)/2
                  Call dDaFile(Lu_A(iIrrep),2,Am,nAm,iAddr_)
               End If
*              Point to Q_k in Qm
               iOff = (kCol-1)*kCol/2
               ip_Q_l=ip_Qm+iOff
            Else If (kCol.gt.nMem) Then
*              Use special scratch for A_k
               ip_A_l=ip_A_k
               Call dDaFile(Lu_A(iIrrep),2,Work(ip_A_l),kCol,iAddr_)
*              Use special scratch for Q_k
               ip_Q_l=ip_Q_k
            End If
*
            LinDep=2
            Call Inv_Cho_Factor(Work(ip_A_l),kCol,
     &                          Am,Qm,nMem,
     &                          Lu_A(iIrrep),Lu_Q(iIrrep),
     &                          Scr,lScr,
     &                          Z,X,ThrQ,
     &                          Work(ip_Q_l),LinDep)

            If (LinDep.ne.0) Then
               Call WarningMessage(2,'Error in Post_2Center_RI')
               Write(6,*) 'Inv_Cho_Factor found linear dependence!'
               Call Abend()
            End If
*
*           Write the new A/Q-vector to file
*
            iAddr_=iAddr
            If (kCol.eq.nMem) Then
               nQm=kCol*(kCol+1)/2
               Call dDaFile(Lu_Q(iIrrep),1,Qm,nQm,iAddr )
               Call dDaFile(Lu_A(iIrrep),1,Am,nQm,iAddr_)
            Else If (kCol.gt.nMem) Then
               Call dDaFile(Lu_Q(iIrrep),1,Work(ip_Q_l),kCol,iAddr )
               Call dDaFile(Lu_A(iIrrep),1,Work(ip_A_l),kCol,iAddr_)
            End If
*
         End Do
*
         Q_l=>Null()
         A_l=>Null()
         Call mma_deallocate(Q_k)
         Call mma_deallocate(A_k)
         Call mma_deallocate(Qm)
         Call mma_deallocate(Am)
         Call mma_deallocate(X)
         Call mma_deallocate(Z)
         Call mma_deallocate(Scr)
         Call DaClos(Lu_A(iIrrep))
 777     Continue
      End Do ! iIrrep
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Sort the Q-matrix back to the original order.
*
      Call GetMem('MemMax','Max','Real',iDummy,MaxMem2)
*
      If (Force_Out_of_Core) MaxMem2=2*nBfnTot
      lScr=Min(MaxMem2,Max(nBfn2,2*nBfnTot))
      Call mma_allocate(Scr,lScr,Label='Scr')
*
      Call SORT_mat(irc,Work(ipA_Diag),nDmA,nDmB,iDiag,nIrrep,
     &                  LU_Q,'Restore',lScr,Scr)
*
*     Note: after the 'Restore' call to Sort_mat, the Q-matrix is
*           no longer stored as upper-triangular but as RECTANGULAR
*           (nDmA,nDmB). The column index is still pivoted.
*
      Call mma_deallocate(Scr)
      Call mma_deallocate(iDiag)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Return
      End
