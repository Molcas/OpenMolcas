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
      SubRoutine PipekMezey_Iter(Functional,CMO,Ovlp,Thrs,ThrRot,
     &                           ThrGrad,
     &                           iTab_Ptr,nBas_per_Atom,nBas_Start,
     &                           Name,
     &                           nBasis,nOrb2Loc,nAtoms,nMxIter,
     &                           Maximisation,Converged,Debug,Silent)
C
C     Author: T.B. Pedersen
C
C     Based on the original routines by Y. Carissan.
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Real*8      CMO(nBasis,*), Ovlp(nBasis,*)
      Integer     iTab_Ptr(nAtoms)
      Integer     nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
      Character*(LENIN8) Name(nBasis)
      Logical     Maximisation, Converged, Debug, Silent
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

C     Initialization (iteration 0).
C     -----------------------------

      If (.not.Silent) Call CWTime(C1,W1)
      nIter=0
      lRmat = nOrb2Loc**2
      Call GetMem('Rmat','Allo','Real',ipRmat,lRmat)
      Call GenerateP(Ovlp,CMO,Name,nBasis,nOrb2Loc,nAtoms,iTab_Ptr,
     &               nBas_per_Atom,nBas_Start,Debug)
      Call ComputeFunc(nAtoms,nOrb2Loc,iTab_Ptr,Functional,Debug)
      Call GetGrad_PM(nAtoms,nOrb2Loc,iTab_Ptr,GradNorm,Work(ipRmat),
     &                Debug)
      OldFunctional=Functional
      FirstFunctional=Functional
      Delta=Functional
      If (.not.Silent) Then
         Call CWTime(C2,W2)
         TimC = C2 - C1
         TimW = W2 - W1
         Write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1),1X,F7.2)')
     &   nIter,Functional,Delta,GradNorm,TimC,TimW,Zero
      End If

C     Iterations.
C     -----------

      l_PACol = 2*nOrb2Loc
      Call GetMem('PACol','Allo','Real',ip_PACol,l_PACol)
      Converged = .False.
      Do While (nIter.lt.nMxIter .and. .not.Converged)
         If (.not.Silent) Call CWTime(C1,W1)
         Call RotateOrb(Ovlp,CMO,Work(ip_PACol),
     &                  nBasis,nAtoms,iTab_Ptr,
     &                  Maximisation,nOrb2Loc,Name,nBas_per_Atom,
     &                  nBas_Start,ThrRot,PctSkp,
     &                  Debug)
         Call ComputeFunc(nAtoms,nOrb2Loc,iTab_Ptr,Functional,Debug)
         Call GetGrad_PM(nAtoms,nOrb2Loc,iTab_Ptr,GradNorm,Work(ipRmat),
     &                   Debug)
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
      Call GetMem('PACol','Free','Real',ip_PACol,l_PACol)
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
            Write(6,'(A,I8)')    'Number of localised orbitals  : ',
     &                                                 nOrb2loc
            Write(6,'(A,1P,D20.10)') 'Value of P before localisation: ',
     &                                                 FirstFunctional
            Write(6,'(A,1P,D20.10)') 'Value of P after localisation : ',
     &                                                 Functional
         End If
      End If

      End
