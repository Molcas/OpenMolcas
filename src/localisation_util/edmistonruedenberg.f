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
      SubRoutine EdmistonRuedenberg(Functional,CMO,Thrs,ThrRot,ThrGrad,
     &                              nBas,nOrb2Loc,nFro,
     &                              nSym,nMxIter,
     &                              Maximisation,Converged,Debug,Silent)
C
C     Author: T.B. Pedersen
C
C     Purpose: Edmiston-Ruedenberg localisation of occupied orbitals.
C
      Implicit Real*8 (a-h,o-z)
      Real*8  CMO(*)
      Integer nBas(nSym), nOrb2Loc(nSym), nFro(nSym)
      Logical Maximisation, Converged, Debug, Silent
#include "WrkSpc.fh"

      Character*18 SecNam
      Parameter (SecNam = 'EdmistonRuedenberg')

      Character*80 Txt

C     Symmetry is NOT allowed.
C     ------------------------

      If (nSym .ne. 1) Then
         Call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
      End If

C     Initializations.
C     ----------------

      Functional = -9.9D9

      nBasT     = nBas(1)
      nOrb2LocT = nOrb2Loc(1)
      nFroT     = nFro(1)

      Converged = .False.

      irc = -1
      FracMem = 0.3d0 ! 30 percent of memory used as vector buffer
      Call Cho_X_Init(irc,FracMem)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I6)') 'Cho_X_Init returned',irc
         Call SysAbendMsg(SecNam,'Cholesky initialization error:',Txt)
      End If

C     Localise orbitals.
C     ------------------

      kOffC = nBasT*nFroT + 1
      Call EdmistonRuedenberg_Iter(Functional,CMO(kOffC),Thrs,ThrRot,
     &                             ThrGrad,
     &                             nBasT,nOrb2LocT,nMxIter,
     &                             Maximisation,Converged,Debug,
     &                             Silent)

C     Finalizations.
C     --------------

      irc = -1
      Call Cho_X_Final(irc)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I6)') 'Cho_X_Final returned',irc
         Call SysAbendMsg(SecNam,'Cholesky finalization error:',Txt)
      End If

      End
