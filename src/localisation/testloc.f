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
      SubRoutine TestLoc(irc)
C
C     Author: T.B. Pedersen
C
C     Purpose: test localisation:
C
C              1) CC^T = XX^T ?
C              2) U=C^TSX , UU^T = diag ?
C              3) If C^TSC=1, X^TSX=1 ?
C
C              where C are the original MOs and X the localised ones.
C
C              Return codes: irc=0 (all OK), irc=1 (failure).
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"
#include "debug.fh"
#include "WrkSpc.fh"

      Character*7 SecNam
      Parameter (SecNam = 'TestLoc')

      Character*8  Label
      Character*80 Txt

      Logical Prnt

      external ddot_
      External iPrintLevel

C     Set return code.
C     ----------------

      irc = 0

C     Set tolerance; allow more deviation for PAOs.
C     ---------------------------------------------

      If (LocPAO) Then
         Tol = 1.0d-4
      Else
         Tol = 1.0d-6
      End If

C     Get the AO overlap matrix.
C     --------------------------

      lOvlp = nBas(1)**2
      lOaux = 4 + nBas(1)*(nBas(1)+1)/2
      Do iSym = 2,nSym
         lOvlp = lOvlp + nBas(iSym)**2
         lOaux = lOaux + nBas(iSym)*(nBas(iSym)+1)/2
      End Do
      Call GetMem('TstOvlp','Allo','Real',ipOvlp,lOvlp)
      Call GetMem('TstOaux','Allo','Real',ipOaux,lOaux)
      jrc = -1
      iOpt = 2
      iComp = 1
      iSyLbl = 1
      Label = 'Mltpl  0'
      Call RdOne(jrc,iOpt,Label,iComp,Work(ipOaux),iSyLbl)
      If (jrc .ne. 0) Then
         Write(Txt,'(A,I4)') 'RdOne returned',jrc
         Call SysAbendMsg(SecNam,'I/O error!',Txt)
      End If
      Prnt=Debug.and.iPrintLevel(-1).ge.5
      kTri = ipOaux
      kSqr = ipOvlp
      Do iSym = 1,nSym
         Call Tri2Rec(Work(kTri),Work(kSqr),nBas(iSym),Prnt)
         kTri = kTri + nBas(iSym)*(nBas(iSym)+1)/2
         kSqr = kSqr + nBas(iSym)**2
      End Do
      Call GetMem('TstOaux','Free','Real',ipOaux,lOaux)

C     Memory allocation.
C     ------------------

      lUmat = nOrb2Loc(1)**2
      lScr  = nBas(1)*nOrb2Loc(1)
      Do iSym = 2,nSym
         lUmat = lUmat + nOrb2Loc(iSym)**2
         lScr  = max(lScr,nBas(iSym)*nOrb2Loc(iSym))
      End Do
      lDen = lOvlp
      Call GetMem('DenC','Allo','Real',ipDenC,lDen)
      Call GetMem('DenX','Allo','Real',ipDenX,lDen)
      Call GetMem('Ddff','Allo','Real',ipDdff,lDen)
      Call GetMem('Scratch','Allo','Real',ipScr,lScr)
      Call GetMem('Umat','Allo','Real',ipUmat,lUmat)

C     Test 1) density.
C     ----------------

      nErr = 0

      kDC = ipDenC
      kDX = ipDenX
      kDd = ipDdff
      kC  = ipMOrig
      kX  = ipCMO
      Do iSym = 1,nSym
         kC1 = kC + nBas(iSym)*nFro(iSym)
         Call GetDens_Localisation(Work(kDC),Work(kC1),
     &                             nBas(iSym),nOrb2Loc(iSym))
         kX1 = kX + nBas(iSym)*nFro(iSym)
         Call GetDens_Localisation(Work(kDX),Work(kX1),
     &                             nBas(iSym),nOrb2Loc(iSym))
         nB2 = nBas(iSym)**2
         Call dCopy_(nB2,Work(kDC),1,Work(kDd),1)
         Call dAXPY_(nB2,-1.0d0,Work(kDX),1,Work(kDd),1)
         xNrm = sqrt(dDot_(nB2,Work(kDd),1,Work(kDd),1))
         If (xNrm .gt. Tol) Then
            Write(6,'(A,A,D16.8,A,I2,A)')
     &      SecNam,': ERROR: ||CC^T - XX^T|| = ',xNrm,' (sym.',
     &      iSym,')'
            nErr = nErr + 1
         End If
         kDC = kDC + nB2
         kDX = kDX + nB2
         kDd = kDd + nB2
         kC  = kC  + nB2
         kX  = kX  + nB2
      End Do
      If (nErr .ne. 0) Then
         irc = 1
         Go To 1 ! exit, after de-allocations
      End If

C     Test 2) U^TU = diag.
C     --------------------

      nErr = 0

      kU = ipUmat
      kO = ipOvlp
      kC = ipMOrig
      kX = ipCMO
      Do iSym = 1,nSym
         kC1 = kC + nBas(iSym)*nFro(iSym)
         kX1 = kX + nBas(iSym)*nFro(iSym)
         Call GetUmat_Localisation(Work(kU),Work(kC1),Work(kO),
     &                             Work(kX1),Work(ipScr),lScr,
     &                             nBas(iSym),nOrb2Loc(iSym))
         nTO = max(nOrb2Loc(iSym),1)
         Call DGEMM_('T','N',
     &              nOrb2Loc(iSym),nOrb2Loc(iSym),nOrb2Loc(iSym),
     &              1.0d0,Work(kU),nTO,Work(kU),nTO,
     &              0.0d0,Work(ipScr),nTO)
         xErr = -9.9d9
         ip0  = ipScr - 1
         Do j = 1,nOrb2Loc(iSym)
            kOff = ip0 + nOrb2Loc(iSym)*(j-1)
            Do i = j+1,nOrb2Loc(iSym)
               Tst  = abs(Work(kOff+i))
               xErr = max(xErr,Tst)
            End Do
         End Do
         If (xErr .gt. Tol) Then
            Write(6,'(A,A,D16.8,A,I2,A)')
     &      SecNam,': ERROR: max. U^TU off-diag. = ',xErr,' (sym.',
     &      iSym,')'
            nErr = nErr + 1
         End If
         nB2 = nBas(iSym)**2
         kU  = kU + nOrb2Loc(iSym)**2
         kO  = kO + nB2
         kC  = kC + nB2
         kX  = kX + nB2
      End Do
      If (nErr .ne. 0) Then
         irc = 1
         Go To 1 ! exit, after de-allocations
      End If

C     Test 3) X^TSX=1.
C     ----------------

      nErr = 0

      kU = ipUmat
      kO = ipOvlp
      kC = ipMOrig
      kX = ipCMO
      Do iSym = 1,nSym
         kC1 = kC + nBas(iSym)*nFro(iSym)
         Call GetUmat_Localisation(Work(kU),Work(kC1),Work(kO),
     &                             Work(kC1),Work(ipScr),lScr,
     &                             nBas(iSym),nOrb2Loc(iSym))
         mErr = 0
         xErr = -9.9d9
         ip0  = kU - 1
         Do j = 1,nOrb2Loc(iSym)
            kOff = ip0 + nOrb2Loc(iSym)*(j-1)
            Tst  = abs(Work(kOff+j)-1.0d0)
            If (Tst .gt. Tol) Then
               mErr = mErr + 1
            End If
            Do i = j+1,nOrb2Loc(iSym)
               Tst  = abs(Work(kOff+i))
               xErr = max(xErr,Tst)
            End Do
         End Do
         If (xErr .gt. Tol) Then
            mErr = mErr + 1
         End If
         If (mErr .eq. 0) Then
            kX1 = kX + nBas(iSym)*nFro(iSym)
            Call GetUmat_Localisation(Work(kU),Work(kX1),Work(kO),
     &                                Work(kX1),Work(ipScr),lScr,
     &                                nBas(iSym),nOrb2Loc(iSym))
            xErr = -9.9d9
            ip0  = kU - 1
            Do j = 1,nOrb2Loc(iSym)
               kOff = ip0 + nOrb2Loc(iSym)*(j-1)
               Tst  = abs(Work(kOff+j)-1.0d0)
               If (Tst .gt. Tol) Then
                  mErr = mErr + 1
               End If
               Do i = j+1,nOrb2Loc(iSym)
                  Tst  = abs(Work(kOff+i))
                  xErr = max(xErr,Tst)
               End Do
            End Do
            If (xErr .gt. Tol) Then
               Write(6,'(A,A,D16.8,A,I2,A)')
     &         SecNam,': ERROR: max. X^TSX off-diag. = ',xErr,' (sym.',
     &         iSym,')'
               mErr = mErr + 1
            End If
            If (mErr .ne. 0) nErr = nErr + 1
         Else
            Write(6,*) SecNam,': original orbitals not orthonormal!'
            Write(6,*) 'Orthonormality test is skipped (sym. ',iSym,')'
         End If
         nB2 = nBas(iSym)**2
         kU  = kU + nOrb2Loc(iSym)**2
         kO  = kO + nB2
         kC  = kC + nB2
         kX  = kX + nB2
      End Do
      If (nErr .ne. 0) Then
         irc = 1
         Go To 1 ! exit, after de-allocations
      End If

C     De-allocations and return.
C     --------------------------

    1 Continue
      Call GetMem('Umat','Free','Real',ipUmat,lUmat)
      Call GetMem('Scratch','Free','Real',ipScr,lScr)
      Call GetMem('Ddff','Free','Real',ipDdff,lDen)
      Call GetMem('DenX','Free','Real',ipDenX,lDen)
      Call GetMem('DenC','Free','Real',ipDenC,lDen)
      Call GetMem('TstOvlp','Free','Real',ipOvlp,lOvlp)

      End
