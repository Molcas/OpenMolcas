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
      SubRoutine Chk_Input(irc)
C
C     Author: T.B. Pedersen
C     Purpose: Check input.
C
C     Return codes:
C       irc = 0: all OK
C       irc < 0: all OK, but nothing to do
C       irc > 0: input error
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"

      Character*9 SecNam
      Parameter (SecNam = 'Chk_Input')

      Logical DoCholesky

      irc = 0
      doCholesky = .False.

      nOrb2LocT = 0
      Do iSym = 1,nSym
         n = nFro(iSym) + nOrb2Loc(iSym)
         If (n.lt.0 .or. n.gt.nOrb(iSym)) Then
            irc = irc + 1
            Write(6,*) SecNam,': nFro + nOrb2Loc out of bounds:'
            Write(6,*) '    iSym     = ',iSym
            Write(6,*) '    nFro     = ',nFro(iSym)
            Write(6,*) '    nOrb2Loc = ',nOrb2Loc(iSym)
            Write(6,*) '    nOrb     = ',nOrb(iSym)
         End If
         If (n .gt. nBas(iSym)) Then
            irc = irc + 1
            Write(6,*) SecNam,': nFro + nOrb2Loc > nBas:'
            Write(6,*) '    iSym     = ',iSym
            Write(6,*) '    nFro     = ',nFro(iSym)
            Write(6,*) '    nOrb2Loc = ',nOrb2Loc(iSym)
            Write(6,*) '    nBas     = ',nBas(iSym)
         End If
         nOrb2LocT = nOrb2LocT + nOrb2Loc(iSym)
      End Do
      If (nOrb2LocT .eq. 0) Then
         irc = -1
         Return
      End If

      If (LocModel.lt.0 .or. LocModel.gt.nLocModel) Then
         Write(6,*) SecNam,': LocModel must satisfy 0 <= LocModel <= ',
     &              nLocModel
         Write(6,*) '    LocModel = ',LocModel
         irc = irc + 1
      End If

      If (LocModel .eq. 4) Then
         Call DecideOnCholesky(doCholesky)
         If (.not.doCholesky) Then
            Call SysAbendMsg(SecNam,
     &           'Edmiston-Ruedenberg localisation not possible:',
     &           'Cholesky integrals required!')
         End If
      End If

      If (EvalER) Then
         Call DecideOnCholesky(doCholesky)
         If (.not.doCholesky) Then
            Write(6,*) SecNam,': evaluation of ER functional requires',
     &                 ' Cholesky decomposition of ERIs!'
            Write(6,*) 'Evaluation of ER functional is cancelled...'
            EvalER = .False.
         End If
      End If

      If (Analysis .and. .not.Test_Localisation) Then
         Test_Localisation = .True.
      End If

      End
