!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SubRoutine Cho_SetShP2Q_2(irc,iLoc,iShlAB,nAB)
!
!     Purpose: set mapping from shell pair iShlAB to qualified
!              columns within current reduced set (stored at location
!              iLoc = 2 or 3).
!              If a non-zero code (irc) is returned, nothing has been
!              set!!
!
      use ChoArr, only: iSP2F, nBstSh, iShP2Q
      use ChoSwp, only: iQuAB, IndRSh, IndRed
      Implicit Real*8 (a-h,o-z)
      Integer nAB(*)
#include "cholesky.fh"

!     Check allocations.
!     ------------------

      Call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.True.)
      If (iShlA .eq. iShlB) Then
         NumAB = nBstSh(iShlA)*(nBstSh(iShlA)+1)/2
      Else
         NumAB = nBstSh(iShlA)*nBstSh(iShlB)
      End If
      lTst = 2*NumAB
      l_iShP2Q = 0
      If (Allocated(iShP2Q)) l_iShP2Q=SIZE(iShP2Q)
      If (l_iShP2Q.lt.1 .or. l_iShP2Q.lt.lTst) Then
         irc = 102
         Return
      End If

!     Check iLoc.
!     -----------

      If (iLoc.lt.2 .or. iLoc.gt.3) Then
         irc = 104
         Return
      End If

!     Set mapping array.
!     iShP2Q(1,AB) = index among qualified, symmetry reduced.
!     iShP2Q(2,AB) = symmetry block.
!     Zeros are returned if the element AB is not qualified.
!     -------------------------------------------------------

      iShP2Q(:,1:NumAB)=0
      Call iZero(nAB,nSym)

      Do iSym = 1,nSym
         Do iQ = 1,nQual(iSym)
            iAB = iQuAB(iQ,iSym) ! addr in curr. red. set
            jAB = IndRed(iAB,iLoc)  ! addr in 1st red. set
            kShlAB = IndRSh(jAB)  ! shell pair (full)
            If (kShlAB .eq. iSP2F(iShlAB)) Then
               kAB = IndRed(jAB,1) ! addr in full shell pair
               iShP2Q(1,kAB)   = iQ
               iShP2Q(2,kAB) = iSym
               nAB(iSym) = nAB(iSym) + 1
            End If
         End Do
      End Do

!     Set return code 0: all ok!
!     --------------------------

      irc = 0

      End
