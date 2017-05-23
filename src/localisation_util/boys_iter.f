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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine Boys_Iter(Functional,CMO,Thrs,ThrRot,ThrGrad,
     &                     ipLbl_AO,ipLbl,
     &                     nBas,nOrb2Loc,nComp,nMxIter,
     &                     Maximisation,Converged,Debug,Silent)
C
C     Author: T.B. Pedersen
C
C     Purpose: Boys localisation of orbitals.
C
      Implicit Real*8 (a-h,o-z)
      Real*8  CMO
      Integer ipLbl_AO(nComp), ipLbl(nComp)
      Integer nBas, nOrb2Loc
      Logical Maximisation, Converged, Debug, Silent
#include "real.fh"
#include "WrkSpc.fh"

C     Print iteration table header.
C     -----------------------------

      If (.not.Silent) Then
         Write(6,'(//,1X,A,A,/,1X,A,A)')
     &   '                                                        CPU ',
     &   '      Wall',
     &   'nIter       Functional P        Delta     Gradient     (sec)',
     &   '     (sec) %Screen'
      End If

C     Initialization (iter 0).
C     ------------------------

      If (.not.Silent) Call CWTime(C1,W1)
      nIter=0
      Converged = .False.
      lRmat = nOrb2Loc**2
      Call GetMem('Rmat','Allo','Real',ipRmat,lRmat)
      Call GenerateB(CMO,nBas,nOrb2Loc,ipLbl_AO,ipLbl,nComp,Debug)
      Call ComputeFuncB2(nOrb2Loc,ipLbl,nComp,Functional,Debug)
      Call GetGrad_Boys(nOrb2Loc,ipLbl,nComp,Work(ipRmat),GradNorm,
     &                  Debug)
      OldFunctional=Functional
      FirstFunctional=Functional
      Delta=Functional
      If (.not. Silent) Then
        Call CWTime(C2,W2)
        TimC = C2 - C1
        TimW = W2 - W1
        Write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1),1X,F7.2)')
     &  nIter,Functional,Delta,GradNorm,TimC,TimW,Zero
      End If

C     Iterations.
C     -----------

      lCol = 2*nOrb2Loc
      Call GetMem('Col','Allo','Real',ipCol,lCol)
      Do While (nIter.lt.nMxIter .and. .not.Converged)
         If (.not.Silent) Call CWTime(C1,W1)
         Call RotateOrbB(CMO,Work(ipCol),ipLbl,nComp,nBas,nOrb2Loc,
     &                   Maximisation,ThrRot,PctSkp,Debug)
         Call ComputeFuncB2(nOrb2Loc,ipLbl,nComp,Functional,Debug)
         Call GetGrad_Boys(nOrb2Loc,ipLbl,nComp,Work(ipRmat),GradNorm,
     &                     Debug)
         nIter=nIter+1
         Delta=Functional-OldFunctional
         OldFunctional=Functional
         If (.not.Silent) Then
            Call CWTime(C2,W2)
            TimC = C2 - C1
            TimW = W2 - W1
            Write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1),1X,F7.2)')
     &      nIter,Functional,Delta,GradNorm,TimC,TimW,PctSkp
         End If
         Converged=GradNorm.le.ThrGrad .and. abs(Delta).le.Thrs
      End Do
      Call GetMem('Col','Free','Real',ipCol,lCol)
      Call GetMem('Rmat','Free','Real',ipRmat,lRmat)

C     Print convergence message.
C     --------------------------

      If (.not.Silent) Then
         If (.not.Converged) Then
            Write(6,'(/,A,I4,A)')
     &      'No convergence after',nIter,' iterations.'
         Else
            Write(6,'(/,A,I4,A)')
     &      'Convergence after',nIter,' iterations.'
            Write(6,*)
            Write(6,'(A,1X,I4)')    'Number of localised orbitals  :',
     &                              nOrb2Loc
            Write(6,'(A,1X,1P,D20.10)')
     &                              'Value of P before localisation:',
     &                              FirstFunctional
            Write(6,'(A,1X,1P,D20.10)')
     &                              'Value of P after localisation :',
     &                              Functional
         End If
      End If

      End
