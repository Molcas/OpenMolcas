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
     &                           ThrGrad,PA,nBas_per_Atom,nBas_Start,
     &                           Name,nBasis,nOrb2Loc,nAtoms,nMxIter,
     &                           Maximisation,Converged,Debug,Silent)
C
C     Author: T.B. Pedersen
C
C     Based on the original routines by Y. Carissan.
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8      CMO(nBasis,*), Ovlp(nBasis,*)
      Real*8      PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Integer     nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
      Character*(LENIN8) Name(nBasis)
      Logical     Maximisation, Converged, Debug, Silent
      Real*8, Allocatable:: RMat(:,:), PACol(:,:)

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
      Call mma_Allocate(RMat,nOrb2Loc,nOrb2Loc,Label='RMat')
      Call GenerateP(Ovlp,CMO,Name,nBasis,nOrb2Loc,nAtoms,
     &               nBas_per_Atom,nBas_Start,PA,Debug)
      Call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
      Call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug)
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

      Call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')
      Converged = .False.
      Do While (nIter.lt.nMxIter .and. .not.Converged)
         If (.not.Silent) Call CWTime(C1,W1)
         Call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,
     &                  Maximisation,nOrb2Loc,Name,nBas_per_Atom,
     &                  nBas_Start,ThrRot,PctSkp,Debug)
         Call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
         Call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug)
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
      Call mma_Deallocate(PACol)
      Call mma_Deallocate(RMat)

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
