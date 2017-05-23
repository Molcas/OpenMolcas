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
      SubRoutine OrthoPAO_Localisation(X,nBas,nFro,nOrb2Loc,nSym,nPass,
     &                                 Test)
C
C     Thomas Bondo Pedersen, December 2005.
C     - revised January 2006.
C
C     Purpose: orthonormalization of Cholesky PAOs according to
C
C              V = X^T*S*X
C              X <- X*V^(-1/2)
C
C              where S is the AO overlap matrix.
C              The orthonormalization is carried out nPass times.
C              After this routine, X will satisfy X^T*S*X=1.
C
C     NOTE: X is assumed to contain all orbitals!!
C
      Implicit Real*8 (a-h,o-z)
      Real*8  X(*)
      Integer nBas(nSym), nFro(nSym), nOrb2Loc(nSym)
      Logical Test
#include "WrkSpc.fh"

      Character*21 SecNam
      Parameter (SecNam = 'OrthoPAO_Localisation')

      Parameter (Tol = 1.0d-10)

      external ddot_

C     Check for quick return.
C     -----------------------

      If (nPass .lt. 1) Return

C     Read S from disk, stored as full square.
C     ----------------------------------------

      l_S = nBas(1)**2
      Do iSym = 2,nSym
         l_S = l_S + nBas(iSym)**2
      End Do
      Call GetMem('S','Allo','Real',ip_S,l_S)
      Call GetOvlp_Localisation(Work(ip_S),'Sqr',nBas,nSym)

C     Allocations.
C     ------------

      nBasMax = nBas(1)
      nO2LMax = nOrb2Loc(1)
      Do iSym = 2,nSym
         nBasMax = max(nBasMax,nBas(iSym))
         nO2LMax = max(nO2LMax,nOrb2Loc(iSym))
      End Do

      l_V = nO2LMax**2
      l_VSqrt = l_V
      l_VISqrt = l_V
      l_Scr = 2*(nBasMax**2) + nBasMax*(nBasMax+1)/2  ! needed in SqrtMt
      Call GetMem('V','Allo','Real',ip_V,l_V)
      Call GetMem('VSqrt','Allo','Real',ip_VSqrt,l_VSqrt)
      Call GetMem('VISqrt','Allo','Real',ip_VISqrt,l_VISqrt)
      Call GetMem('Scr','Allo','Real',ip_Scr,l_Scr)

C     Orthonormalization passes.
C     --------------------------

      Do iPass = 1,nPass
         kX = 1
         kS = ip_S
         Do iSym = 1,nSym

C           Set pointer to PAO part of X.
C           ------------------------------

            kOffX = kX + nBas(iSym)*nFro(iSym)

C           Compute V = X^T*S*X.
C           --------------------

            Call GetUmat_Localisation(Work(ip_V),
     &                                X(kOffX),Work(kS),X(kOffX),
     &                                Work(ip_Scr),l_Scr,
     &                                nBas(iSym),nOrb2Loc(iSym))

C           Compute V^(-1/2).
C           -----------------

            iTask = 2 ! compute sqrt as well as inverse sqrt
            Call SqrtMt(Work(ip_V),nOrb2Loc(iSym),iTask,
     &                  Work(ip_VSqrt),Work(ip_VISqrt),Work(ip_Scr))

C           Compute orthonormal X <- X*V^(-1/2).
C           ------------------------------------

            nB = max(nBas(iSym),1)
            nO2L = max(nOrb2Loc(iSym),1)
            Call dCopy_(nBas(iSym)*nOrb2Loc(iSym),X(kOffX),1,
     &                                           Work(ip_Scr),1)
         Call DGEMM_('N','N',nBas(iSym),nOrb2Loc(iSym),nOrb2Loc(iSym),
     &                 1.0d0,Work(ip_Scr),nB,Work(ip_VISqrt),nO2L,
     &                 0.0d0,X(kOffX),nB)

C           Update pointers.
C           ----------------

            kX = kX + nBas(iSym)**2
            kS = kS + nBas(iSym)**2

         End Do
      End Do

C     Test orthonormalization (i.e. V=1?).
C     ------------------------------------

      If (Test) Then
         kX = 1
         kS = ip_S
         nErr = 0
         Do iSym = 1,nSym
            kOffX = kX + nBas(iSym)*nFro(iSym)
            Call GetUmat_Localisation(Work(ip_V),
     &                                X(kOffX),Work(kS),X(kOffX),
     &                                Work(ip_Scr),l_Scr,
     &                                nBas(iSym),nOrb2Loc(iSym))
            kOff = ip_V - 1
            Do i = 1,nOrb2Loc(iSym)
               Work(kOff+nOrb2Loc(iSym)*(i-1)+i) =
     &              Work(kOff+nOrb2Loc(iSym)*(i-1)+i) - 1.0d0
            End Do
            xNrm = sqrt(dDot_(nOrb2Loc(iSym)**2,Work(ip_V),1,
     &                                         Work(ip_V),1))
            If (xNrm .gt. Tol) Then
               Write(6,'(A,A,D16.8,A,I2,A)')
     &         SecNam,': ERROR: ||X^TSX - 1|| = ',xNrm,' (sym.',
     &         iSym,')'
               nErr = nErr + 1
            End If
            kX = kX + nBas(iSym)**2
            kS = kS + nBas(iSym)**2
         End Do
         If (nErr .ne. 0) Then
            Write(6,*) SecNam,': failure after ',nPass,' passes'
            Call SysAbendMsg(SecNam,'Orthonormalization failure!',' ')
         End If
      End If

C     De-allocations.
C     ---------------

      Call GetMem('Scr','Free','Real',ip_Scr,l_Scr)
      Call GetMem('VISqrt','Free','Real',ip_VISqrt,l_VISqrt)
      Call GetMem('VSqrt','Free','Real',ip_VSqrt,l_VSqrt)
      Call GetMem('V','Free','Real',ip_V,l_V)
      Call GetMem('S','Free','Real',ip_S,l_S)

      End
