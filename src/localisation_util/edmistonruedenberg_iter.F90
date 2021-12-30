!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************
      SubRoutine EdmistonRuedenberg_Iter(Functional,CMO,Thrs,ThrRot,    &
     &                                   ThrGrad,                       &
     &                                   nBasis,nOrb2Loc,nMxIter,       &
     &                                   Maximisation,Converged,Debug,  &
     &                                   Silent)
!
!     Thomas Bondo Pedersen, November 2005.
!
!     Purpose: ER localisation of orbitals.
!
!     The optimization algorithm is the "generalized eta step" of
!     Subotnik, Shao, Liang, and Head-Gordon, JCP 121, 9220 (2004).
!
!     Redundant arguments that might be used at a later stage:
!       ThrRot [might be used for DIIS]
!       Maximisation [might be used for Jacobi sweeps]
!
!     Note that two-electron integrals (Cholesky decomposed) must be
!     available and appropriately set up when calling this routine.
!
      Implicit Real*8 (a-h,o-z)
      Real*8  CMO(nBasis,nOrb2Loc)
      Logical Maximisation, Converged, Debug, Silent
#include "WrkSpc.fh"

      Character*23 SecNam
      Parameter (SecNam = 'EdmistonRuedenberg_Iter')

      Logical Timing

      If (Debug) Then
         Write(6,*) SecNam,'[debug]: Maximisation: ',Maximisation
         Write(6,*) SecNam,'[debug]: ThrRot      : ',ThrRot
      End If

!     Print iteration table header.
!     -----------------------------

      If (.not.Silent) Then
         Write(6,'(//,1X,A,A,/,1X,A,A)')                                &
     &   '                                                        CPU ',&
     &   '      Wall',                                                  &
     &   'nIter      Functional ER        Delta     Gradient     (sec)',&
     &   '     (sec)'
      End If

!     Initialization.
!     ---------------

      Converged = .False.
      Timing = Debug

      lRmat = nOrb2Loc**2
      Call GetMem('Rmat','Allo','Real',ipRmat,lRmat)

!     Iteration 0.
!     ------------

      If (.not.Silent) Call CWTime(C1,W1)
      nIter = 0
      Functional = 0.0d0
      Call GetGrad_ER(Functional,GradNorm,Work(ipRmat),CMO,             &
     &                nBasis,nOrb2Loc,Timing)
      OldFunctional = Functional
      FirstFunctional = Functional
      Delta = Functional
      If (.not.Silent) Then
         Call CWTime(C2,W2)
         TimC = C2 - C1
         TimW = W2 - W1
         Write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1))')             &
     &   nIter,Functional,Delta,GradNorm,TimC,TimW
      End If

!     Iterations.
!     -----------

      Do While (nIter.lt.nMxIter .and. .not.Converged)
         If (.not.Silent) Call CWTime(C1,W1)
         Call RotateOrb_ER(Work(ipRmat),CMO,nBasis,nOrb2Loc,Debug)
         Call GetGrad_ER(Functional,GradNorm,Work(ipRmat),              &
     &                   CMO,nBasis,nOrb2Loc,Timing)
         nIter = nIter + 1
         Delta = Functional - OldFunctional
         OldFunctional = Functional
         If (.not.Silent) Then
            Call CWTime(C2,W2)
            TimC = C2 - C1
            TimW = W2 - W1
            Write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1))')          &
     &      nIter,Functional,Delta,GradNorm,TimC,TimW
         End If
         Converged=GradNorm.le.ThrGrad .and. abs(Delta).le.Thrs
      End Do

!     Print convergence message.
!     --------------------------

      If (.not.Silent) Then
         If (.not.Converged) Then
            Write(6,'(/,A,I4,A)')                                       &
     &      'No convergence after',nIter,' iterations.'
         Else
            Write(6,'(/,A,I4,A)')                                       &
     &      'Convergence after',nIter,' iterations.'
            Write(6,*)
            Write(6,'(A,I8)')    'Number of localised orbitals  : ',    &
     &                                                 nOrb2loc
            Write(6,'(A,F12.8)') 'Value of P before localisation: ',    &
     &                                                 FirstFunctional
            Write(6,'(A,F12.8)') 'Value of P after localisation : ',    &
     &                                                 Functional
         End If
      End If

!     Finalization.
!     -------------

      Call GetMem('Rmat','Free','Real',ipRmat,lRmat)

      End
