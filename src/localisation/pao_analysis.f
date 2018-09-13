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
* Copyright (C) 2005,2006, Thomas Bondo Pedersen                       *
************************************************************************
      SubRoutine PAO_Analysis(D,R,X)
C
C     Thomas Bondo Pedersen, December 2005.
C     - revised January 2006 (Thomas Bondo Pedersen).
C
C     Purpose: test and analysis of Cholesky PAOs.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 D(*), R(*), X(*)
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"

      l_S = nBas(1)**2
      Do iSym = 2,nSym
         l_S = l_S + nBas(iSym)**2
      End Do
      Call GetMem('S','Allo','Real',ip_S,l_S)
      Call GetOvlp_Localisation(Work(ip_S),'Sqr',nBas,nSym)
      Call PAO_Ana1(D,R,X,Work(ipMOrig),Work(ip_S),
     &              nBas,nFro,nOrb2Loc,nSym,Name,nAtoms,AnaNrm)
      Call GetMem('S','Free','Real',ip_S,l_S)

      End
      SubRoutine PAO_Ana1(D,R,X,C,S,nBas,nFro,nOrb2Loc,nSym,
     &                    Nam,nAtoms,AnaNrm)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Integer nBas(nSym), nFro(nSym), nOrb2Loc(nSym)
      Real*8  D(*), R(*), X(*), C(*), S(*)
      Character*(LENIN8) Nam(*)
      Character*3 AnaNrm
#include "WrkSpc.fh"

      Character*8 SecNam
      Parameter (SecNam = 'PAO_Ana1')

      Character*14 FilNam
      Logical Debug

      external ddot_

      Parameter (Tol = 1.0d-10)

C     Initialization.
C     ---------------

      Debug = .False.

C     Tests:
C     0) RR^T = XX^T = D (check of decomposition)
C     1) Co^TSX = 0 (check that X is orthogonal to the orthogonal
C                    complement space).
C     2) C^TSX is non-singular (check that X spans primary space)
C     ===========================================================

      Write(6,*) 'Testing PAOs...'
      nErr = 0

C     Test 0.
C     -------

      l_Tst = nBas(1)**2
      Do iSym = 2,nSym
         l_Tst = max(l_Tst,nBas(iSym)**2)
      End Do
      Call GetMem('TstDen','Allo','Real',ip_Tst,l_Tst)

      kOff = 1
      Do iSym = 1,nSym
         nB = nBas(iSym)
         nB2 = nB**2
         nF = nFro(iSym)
         nO = nOrb2Loc(iSym)
         Call GetDens_Localisation(Work(ip_Tst),R(kOff),nB,nB)
         Call dAXPY_(nB2,-1.0d0,D(kOff),1,Work(ip_Tst),1)
         If (nB2 .gt. 0) Then
            xB2 = 1.0d0/dble(nB2)
            xRMS = sqrt(xB2*dDot_(nB2,Work(ip_Tst),1,Work(ip_Tst),1))
         Else
            xRMS = 0.0d0
         End If
         If (xRMS .gt. Tol) Then
            Write(6,'(A,A,D16.8,A,I2,A)')
     &      SecNam,': ERROR: RMS(D-RR^T) = ',xRMS,' (sym.',
     &      iSym,')'
            nErr = nErr + 1
         End If
         kOffX = kOff + nB*nF
         Call GetDens_Localisation(Work(ip_Tst),X(kOffX),nB,nO)
         Call dAXPY_(nB2,-1.0d0,D(kOff),1,Work(ip_Tst),1)
         If (nB2 .gt. 0) Then
            xB2 = 1.0d0/dble(nB2)
            xRMS = sqrt(xB2*dDot_(nB2,Work(ip_Tst),1,Work(ip_Tst),1))
         Else
            xRMS = 0.0d0
         End If
         If (xRMS .gt. Tol) Then
            Write(6,'(A,A,D16.8,A,I2,A)')
     &      SecNam,': ERROR: RMS(D-XX^T) = ',xRMS,' (sym.',
     &      iSym,')'
            nErr = nErr + 1
         End If
         kOff = kOff + nB2
      End Do

C     Compute SX.
C     -----------

      l_SX = nBas(1)*nOrb2Loc(1)
      Do iSym = 2,nSym
         l_SX = l_SX + nBas(iSym)*nOrb2Loc(iSym)
      End Do
      Call GetMem('SX','Allo','Real',ip_SX,l_SX)

      kOff = 1
      kSX  = ip_SX
      Do iSym = 1,nSym
         nB = max(nBas(iSym),1)
         kX = kOff + nBas(iSym)*nFro(iSym)
         Call DGEMM_('N','N',nBas(iSym),nOrb2Loc(iSym),nBas(iSym),
     &              1.0d0,S(kOff),nB,X(kX),nB,
     &              0.0d0,Work(kSX),nB)
         kSX  = kSX  + nBas(iSym)*nOrb2Loc(iSym)
         kOff = kOff + nBas(iSym)*nBas(iSym)
      End Do

C     Test 1.
C     First frozen part of orthogonal complement, then the rest.
C     ----------------------------------------------------------

      kC   = 1
      kSX  = ip_SX
      kTst = ip_Tst
      Do iSym = 1,nSym
         nB = max(nBas(iSym),1)
         nF = max(nFro(iSym),1)
         nRest = nBas(iSym) - nFro(iSym) - nOrb2Loc(iSym)
         nR = max(nRest,1)
         Call DGEMM_('T','N',nFro(iSym),nOrb2Loc(iSym),nBas(iSym),
     &              1.0d0,C(kC),nB,Work(kSX),nB,
     &              0.0d0,Work(kTst),nF)
         nFO = nFro(iSym)*nOrb2Loc(iSym)
         If (nFO .gt. 0) Then
            xFO = 1.0d0/dble(nFO)
            xRMS = sqrt(xFO*dDot_(nFO,Work(kTst),1,Work(kTst),1))
         Else
            xRMS = 0.0d0
         End If
         If (xRMS .gt. Tol) Then
            Write(6,'(A,A,D16.8,A,I2,A)')
     &      SecNam,': ERROR: RMS(Co^TSX [Frozen]) = ',xRMS,' (sym.',
     &      iSym,')'
            nErr = nErr + 1
         End If
         kCR = kC + nBas(iSym)*(nFro(iSym)+nOrb2Loc(iSym))
         Call DGEMM_('T','N',nRest,nOrb2Loc(iSym),nBas(iSym),
     &              1.0d0,C(kCR),nB,Work(kSX),nB,
     &              0.0d0,Work(kTst),nR)
         nRO = nRest*nOrb2Loc(iSym)
         If (nRO .gt. 0) Then
            xRO = 1.0d0/dble(nRO)
            xRMS = sqrt(xRO*dDot_(nRO,Work(kTst),1,Work(kTst),1))
         Else
            xRMS = 0.0d0
         End If
         If (xRMS .gt. Tol) Then
            Write(6,'(A,A,D16.8,A,I2,A)')
     &      SecNam,': ERROR: RMS(Co^TSX [Rest]) = ',xRMS,' (sym.',
     &      iSym,')'
            nErr = nErr + 1
         End If
         kC  = kC  + nBas(iSym)*nBas(iSym)
         kSX = kSX + nBas(iSym)*nOrb2Loc(iSym)
      End Do

C     Test 2.
C     -------

      kC   = 1
      kSX  = ip_SX
      kTst = ip_Tst
      Do iSym = 1,nSym
         nB = max(nBas(iSym),1)
         nO = max(nOrb2Loc(iSym),1)
         kCO = kC + nBas(iSym)*nFro(iSym)
         Call DGEMM_('T','N',nOrb2Loc(iSym),nOrb2Loc(iSym),nBas(iSym),
     &              1.0d0,C(kCO),nB,Work(kSX),nB,
     &              0.0d0,Work(kTst),nO)
         l_EigR = nOrb2Loc(iSym)
         l_EigI = nOrb2Loc(iSym)
         Call GetMem('EigR','Allo','Real',ip_EigR,l_EigR)
         Call GetMem('EigI','Allo','Real',ip_EigI,l_EigI)
         iGetVecs = 0
         Call Diag_Localisation(Work(kTst),Work(ip_EigR),Work(ip_EigI),
     &                          nOrb2Loc(iSym),iGetVecs)
         kOffR = ip_EigR - 1
         kOffI = ip_EigI - 1
         Do i = 1,nOrb2Loc(iSym)
            xNrm = sqrt(Work(kOffR+i)**2+Work(kOffI+i)**2)
            If (xNrm .lt. Tol) Then
               Write(6,'(A,A,I6,A,A,D16.8,A,I2,A)')
     &         SecNam,': ERROR: ||eigenvalue',i,'|| of ',
     &         'C^TSX = ',xNrm,' (sym.',iSym,')'
               nErr = nErr + 1
            End If
         End Do
         Call GetMem('EigI','Free','Real',ip_EigI,l_EigI)
         Call GetMem('EigR','Free','Real',ip_EigR,l_EigR)
         kC  = kC  + nBas(iSym)*nBas(iSym)
         kSX = kSX + nBas(iSym)*nOrb2Loc(iSym)
      End Do

      Call GetMem('SX','Free','Real',ip_SX,l_SX)
      Call GetMem('TstDen','Free','Real',ip_Tst,l_Tst)

C     Overall check.
C     --------------

      If (nErr .ne. 0) Then
         Call SysAbendMsg(SecNam,'PAO localization failed!',' ')
      Else
         Write(6,*) '...OK!'
      End If

C     Analysis.
C     This part does not work with symmetry.
C     ======================================

      If (nSym .ne. 1) Then
         Write(6,*)
         Write(6,*) SecNam,': symmetry not implemented for analysis',
     &              ' of PAOs. Sorry!'
         Return
      End If

C     Allocate atom based density and CMO matrices.
C     ---------------------------------------------

      lDAt = nAtoms**2
      lRAt = nAtoms*nBas(1)
      lXAt = nAtoms*nOrb2Loc(1)
      Call GetMem('DAt','Allo','Real',ipDAt,lDAt)
      Call GetMem('RAt','Allo','Real',ipRAt,lRAt)
      Call GetMem('XAt','Allo','Real',ipXAt,lXAt)

C     Allocate and get index arrays for basis functions per atom.
C     -----------------------------------------------------------

      l_nBas_per_Atom = nAtoms
      l_nBas_Start    = nAtoms
      Call GetMem('nB_per_Atom','Allo','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Allo','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                 Nam,nBas(1),nAtoms,Debug)

C     Generate bitmaps.
C     -----------------

      Call GetAt_Localisation(D,nBas(1),nBas(1),Work(ipDAt),nAtoms,2,
     &                        iWork(ip_nBas_per_Atom),
     &                        iWork(ip_nBas_Start),AnaNrm)
      Call GetAt_Localisation(R,nBas(1),nBas(1),Work(ipRAt),nAtoms,1,
     &                        iWork(ip_nBas_per_Atom),
     &                        iWork(ip_nBas_Start),AnaNrm)
      kOffX = nBas(1)*nFro(1) + 1
      Call GetAt_Localisation(X(kOffX),nBas(1),nOrb2Loc(1),
     &                        Work(ipXAt),nAtoms,1,
     &                        iWork(ip_nBas_per_Atom),
     &                        iWork(ip_nBas_Start),AnaNrm)
      Write(FilNam,'(A)') 'PAO_Dnsty1.bmp'
      Call GenBMp_Loc(Work(ipDat),nAtoms,nAtoms,FilNam,'r')
      Write(FilNam,'(A)') 'PAO_LnDep1.bmp'
      Call GenBMp_Loc(Work(ipRat),nAtoms,nBas(1),FilNam,'r')
      Write(FilNam,'(A)') 'PAO_Chole1.bmp'
      Call GenBMp_Loc(Work(ipXat),nAtoms,nOrb2Loc(1),FilNam,'r')
      Write(6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

C     Deallocations.
C     --------------

      Call GetMem('nB_Start','Free','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call GetMem('nB_per_Atom','Free','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('XAt','Free','Real',ipXAt,lXAt)
      Call GetMem('RAt','Free','Real',ipRAt,lRAt)
      Call GetMem('DAt','Free','Real',ipDAt,lDAt)

      End
