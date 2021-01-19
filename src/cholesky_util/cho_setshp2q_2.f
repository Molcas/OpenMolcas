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
      SubRoutine Cho_SetShP2Q_2(irc,iLoc,iShlAB,nAB)
C
C     Purpose: set mapping from shell pair iShlAB to qualified
C              columns within current reduced set (stored at location
C              iLoc = 2 or 3).
C              If a non-zero code (irc) is returned, nothing has been
C              set!!
C
      use ChoArr, only: iSP2F, nBstSh
      use ChoSwp, only: iQuAB, IndRSh, IndRed
#include "implicit.fh"
      Integer nAB(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "chosew.fh"
#include "WrkSpc.fh"

C     Check allocations.
C     ------------------

      Call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.True.)
      If (iShlA .eq. iShlB) Then
         NumAB = nBstSh(iShlA)*(nBstSh(iShlA)+1)/2
      Else
         NumAB = nBstSh(iShlA)*nBstSh(iShlB)
      End If
      lTst = 2*NumAB
      If (l_iShP2Q.lt.1 .or. l_iShP2Q.lt.lTst) Then
         irc = 102
         Return
      End If

C     Check iLoc.
C     -----------

      If (iLoc.lt.2 .or. iLoc.gt.3) Then
         irc = 104
         Return
      End If

C     Set mapping array.
C     iShP2Q(1,AB) = index among qualified, symmetry reduced.
C     iShP2Q(2,AB) = symmetry block.
C     Zeros are returned if the element AB is not qualified.
C     -------------------------------------------------------

      Call Cho_iZero(iWork(ip_iShP2Q),lTst)
      Call Cho_iZero(nAB,nSym)

      Do iSym = 1,nSym
         Do iQ = 1,nQual(iSym)
            iAB = iQuAB(iQ,iSym) ! addr in curr. red. set
            jAB = IndRed(iAB,iLoc)  ! addr in 1st red. set
            kShlAB = IndRSh(jAB)  ! shell pair (full)
            If (kShlAB .eq. iSP2F(iShlAB)) Then
               kAB = IndRed(jAB,1) ! addr in full shell pair
               iWork(ip_iShP2Q+2*(kAB-1))   = iQ
               iWork(ip_iShP2Q+2*(kAB-1)+1) = iSym
               nAB(iSym) = nAB(iSym) + 1
            End If
         End Do
      End Do

C     Set return code 0: all ok!
C     --------------------------

      irc = 0

      End
