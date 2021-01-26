************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_P_SetRed_L()
C
C     Purpose: set next local reduced set. The next global reduced set
C              must be available at (global) index array location 2.
C
      use ChoSwp, only: nnBstRSh, nnBstRSh_G
      use ChoSwp, only: iiBstRSh, iiBstRSh_G
      use ChoSwp, only:   IndRed,   IndRed_G
      use ChoArr, only: iL2G, MySP
      Implicit None
#include "cholesky.fh"
#include "choglob.fh"
#include "WrkSpc.fh"

      Character*14 SecNam
      Parameter (SecNam = 'Cho_P_SetRed_L')

      Integer irc, iC, nDim
      Integer i, j, k, i0, k0, l, ll
      Integer iSym, iSP, iShlAB

C     Copy current local reduced set (at location 2) to location 3.
C     -------------------------------------------------------------

      Call Cho_X_RSCopy(irc,2,3)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
         Call Cho_Quit('Error in '//SecNam,104)
      End If

C     Re-initialize local reduced set indices at location 2.
C     ------------------------------------------------------

      nDim = nSym*nnShl
      IndRed(:,2)=0
      Call Cho_iZero(iiBstRSh(:,:,2),nDim)
      Call Cho_iZero(nnBstRSh(:,:,2),nDim)
      Call Cho_iZero(iiBstR(:,2),nSym)
      Call Cho_iZero(nnBstR(:,2),nSym)
      nnBstRT(2) = 0

C     Set local nnBstRSh counter at location 2.
C     -----------------------------------------

      Do iSP = 1,nnShl
         iShlAB = mySP(iSP)
         Do iSym = 1,nSym
            nnBstRSh(iSym,iSP,2) = nnBstRSh_G(iSym,iShlAB,2)
         End Do
      End Do

C     Set remaining reduced set indices (excl. IndRed), location 2.
C     -------------------------------------------------------------

      Call Cho_SetRedInd(2)

C     Set local IndRed to point to local 1st reduced set.
C     ---------------------------------------------------

      iC = 0
      Do iSym = 1,nSym
         Do iSP = 1,nnShl
            iShlAB = mySP(iSP)
            i0 = iiBstR_G(iSym,2) + iiBstRSh_G(iSym,iShlAB,2)
            k0 = iiBstR(iSym,3) + iiBstRSh(iSym,iSP,3)
            Do i = 1,nnBstRSh_G(iSym,iShlAB,2)
               j = IndRed_G(i0+i,2) ! addr in global rs1
               iC = iC + 1
               k = 0
               Do While (k .lt. nnBstRSh(iSym,iSP,3))
                  k = k + 1
                  ll = IndRed(k0+k,3)
                  l = iL2G(ll)
                  If (l .eq. j) Then
                     IndRed(iC,2) = ll
                     k = nnBstRSh(iSym,iSP,3) ! break while loop
                  End If
               End Do
            End Do
         End Do
      End Do

      End
