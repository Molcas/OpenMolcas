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
      SubRoutine Localise_Iterative(irc,Model,Functional)
C
C     Author: T.B. Pedersen
C
C     Purpose: Iterative localisation of orbitals.
C              Models implemented:
C                Pipek-Mezey         [MODEL='PIPE']
C                Boys                [MODEL='BOYS']
C                Edmiston-Ruedenberg [MODEL='EDMI']
C
      Implicit Real*8 (a-h,o-z)
      Character*4 Model
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"
#include "debug.fh"

      Character*18 SecNam
      Parameter (SecNam = 'Localise_Iterative')

      Character*4  myModel
      Character*80 Txt
      Logical Converged

      irc = 0
      Functional = -9.9d9
      Converged  = .False.

C     Generate Cholesky start guess, if requested.
C     --------------------------------------------

      If (ChoStart) Then
         Thrs_Save = Thrs
         Thrs = 1.0d-12
         Call Localise_Noniterative(irc,'Chol',xNrm)
         If (irc .ne. 0) Then
            Write(Txt,'(A,I4)') 'Return code:',irc
            Call SysAbendMsg(SecNam,'Localise_Noniterative failed!',Txt)
         End If
         Thrs = Thrs_Save
      End If

C     Localise.
C     ---------

      myModel = Model
      Call UpCase(myModel)
      If (myModel .eq. 'PIPE') Then
*        If (.not.Silent) Then
            Write(6,'(//,1X,A)') 'Pipek-Mezey localisation'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',Thrs,' (functional)'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',ThrGrad,' (gradient)'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Screening threshold  :',ThrRot,' (orbital rotations)'
            Write(6,'(1X,A,8(1X,I6))')
     &      'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
            Write(6,'(1X,A,8(1X,I6))')
     &      'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
*        End If
         Call PipekMezey(Functional,Work(ipCMO),Thrs,ThrRot,ThrGrad,
     &                   Name,
     &                   nBas,nOrb2Loc,nFro,
     &                   nSym,nAtoms,nMxIter,
     &                   Maximisation,Converged,Debug,Silent)
      Else If (myModel .eq. 'BOYS') Then
*        If (.not.Silent) Then
            Write(6,'(/,1X,A)') 'Boys localisation'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',Thrs,' (functional)'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',ThrGrad,' (gradient)'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Screening threshold  :',ThrRot,' (orbital rotations)'
            Write(6,'(1X,A,8(1X,I6))')
     &      'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
            Write(6,'(1X,A,8(1X,I6))')
     &      'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
*        End If
         Call Boys(Functional,Work(ipCMO),Thrs,ThrRot,ThrGrad,
     &             nBas,nOrb2Loc,nFro,
     &             nSym,nMxIter,
     &             Maximisation,Converged,Debug,Silent)
      Else If (myModel .eq. 'EDMI') Then
*        If (.not.Silent) Then
            Write(6,'(/,1X,A)') 'Edmiston-Ruedenberg localisation'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',Thrs,' (functional)'
            Write(6,'(1X,A,1X,D12.4,A)')
     &      'Convergence threshold:',ThrGrad,' (gradient)'
c           Write(6,'(1X,A,1X,D12.4,A)')
c    &      'Screening threshold  :',ThrRot,' (orbital rotations)'
            Write(6,'(1X,A,8(1X,I6))')
     &      'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
            Write(6,'(1X,A,8(1X,I6))')
     &      'Orbitals to localise :',(nOrb2Loc(iSym),iSym=1,nSym)
*        End If
         Call EdmistonRuedenberg(Functional,Work(ipCMO),Thrs,ThrRot,
     &                           ThrGrad,
     &                           nBas,nOrb2Loc,nFro,
     &                           nSym,nMxIter,
     &                           Maximisation,Converged,Debug,Silent)
      Else
         Write(Txt,'(A,A4)') 'Model = ',Model
         Call SysAbendMsg(SecNam,'Unknown model',Txt)
      End If

      If (.not.Converged) Then
         irc = 1
         Return
      End If

      End
