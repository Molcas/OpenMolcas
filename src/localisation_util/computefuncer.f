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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ComputeFuncER(ERFun,CMO,nBas,nOcc,nFro,nSym,Timing)
C
C     Thomas Bondo Pedersen, November 2005.
C
C     Purpose: compute Edmiston-Ruedenberg functional.
C
C     =====================================================
C        WORKS *ONLY* WITH CHOLESKY DECOMPOSED INTEGRALS
C     =====================================================
C
      Implicit None
      Real*8  ERFun
      Real*8  CMO(*)
      Integer nSym
      Integer nBas(nSym), nOcc(nSym), nFro(nSym)
      Logical Timing
#include "WrkSpc.fh"

      Character*13 SecNam
      Parameter (SecNam = 'ComputeFuncER')

      Character*80 Txt

      Real*8 FracMem

      Integer irc, ipERFun, lERFun, iSym, i, kOff, nFroT
      Integer nOccT(8)

C     Initializations.
C     ----------------

      irc = 0

      FracMem = 0.0d0 ! no buffer allocated
      Call Cho_X_Init(irc,FracMem)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I4)') 'Cho_X_Init returned',irc
         Call SysAbendMsg(SecNam,'Cholesky initialization failure!',Txt)
      End If

C     Check dimensions.
C     -----------------

      Call ERChk_Localisation(irc,nBas,nOcc,nFro,nSym)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I4)') 'ERChk_Localisation returned',irc
         Call SysAbendMsg(SecNam,'Cholesky initialization mismatch!',
     &                    Txt)
      End If

C     Compute ER functional.
C     ----------------------

      nOccT(1) = nOcc(1) + nFro(1)
      Do iSym = 2,nSym
         nOccT(iSym) = nOcc(iSym) + nFro(iSym)
      End Do

      lERFun = nOccT(1)
      nFroT = nFro(1)
      Do iSym = 2,nSym
         lERFun = lERFun + nOccT(iSym)
         nFroT = nFroT + nFro(iSym)
      End Do

      Call GetMem('ERFun','Allo','Real',ipERFun,lERFun)
      ERFun = 0.0d0
      Call EvalERFun(ERFun,Work(ipERFun),CMO,nOccT,nSym,Timing)
      If (nFroT .gt. 0) Then
         kOff = ipERFun - 1
         Do iSym = 1,nSym
            Do i = 1,nFro(iSym)
               ERFun = ERFun - Work(kOff+i)
            End Do
            kOff = kOff + nOccT(iSym)
         End Do
      End If
      Call GetMem('ERFun','Free','Real',ipERFun,lERFun)

C     Finalizations.
C     --------------

      Call Cho_X_Final(irc)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I4)') 'Cho_X_Final returned',irc
         Call SysAbendMsg(SecNam,'Cholesky finalization failure!',Txt)
      End If

      End
      SubRoutine ERChk_Localisation(irc,lnBas,lnOcc,lnFro,lnSym)
      Implicit None
      Integer irc, lnSym
      Integer lnBas(lnSym), lnOcc(lnSym), lnFro(lnSym)
#include "cholesky.fh"
#include "choorb.fh"

      Integer iSym, nTst

      irc = 0

      If (lnSym.lt.1 .or. lnSym.gt.8) Then
         irc = 1
         Return
      End If

      If (lnSym .ne. nSym) Then
         irc = 2
         Return
      End If

      Do iSym = 1,nSym
         If (lnBas(iSym) .ne. nBas(iSym)) Then
            irc = 3
            Return
         End If
         nTst = lnOcc(iSym)+lnFro(iSym)
         If (nTst .gt. nBas(iSym)) Then
            irc = 4
            Return
         End If
      End Do

      End
      SubRoutine EvalERFun(ERFun,ERFunC,CMO,nOcc,nSym,Timing)
      Implicit None
      Real*8  ERFun
      Real*8  ERFunC(*), CMO(*)
      Integer nSym
      Integer nOcc(nSym)
      Logical Timing

      Character*9 SecNam
      Parameter (SecNam = 'EvalERFun')

      Character*80 Txt
      Integer irc

      irc = 0
      Call Cho_Get_ER(irc,CMO,nOcc,ERFunC,ERFun,Timing)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I4)') 'Cho_Get_ER returned',irc
         Call SysAbendMsg(SecNam,'ER evaluation failed!',
     &                    Txt)
      End If

      End
