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
      SubRoutine Cho_SetShP2RS(irc,iLoc,iShlAB,nAB)
!
!     Purpose: set mapping from shell pair iShlAB to current
!              reduced set (stored at location iLoc = 2 or 3).
!              If a non-zero code (irc) is returned, nothing has been
!              set!!
!
      use ChoArr, only: nBstSh, iSP2F, iShP2RS, MySP
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRed
#if defined (_DEBUGPRINT_)
      use ChoSwp, only: IndRSh
#endif
      Implicit Real*8 (a-h,o-z)
      Integer nAB(*)
#include "cholesky.fh"

#if defined (_DEBUGPRINT_)
      Character*13 SecNam
      Parameter (SecNam = 'Cho_SetShP2RS')
#endif

!     Check allocations.
!     ------------------

      Call Cho_InvPck(iSP2F(mySP(iShlAB)),iShlA,iShlB,.True.)
      If (iShlA .eq. iShlB) Then
         NumAB = nBstSh(iShlA)*(nBstSh(iShlA)+1)/2
      Else
         NumAB = nBstSh(iShlA)*nBstSh(iShlB)
      End If
      lTst = 2*NumAB
      l_iShP2RS = 0
      If (Allocated(iShP2RS)) l_iShP2RS = SIZE(iShP2RS)
      If (l_iShP2RS.lt.1 .or. l_iShP2RS.lt.lTst) Then
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
!     iShP2RS(1,AB) = index in current reduced set, symmetry reduced.
!     iShP2RS(2,AB) = symmetry block.
!     Zeros are returned if the element AB is not a member of the
!     current reduced set.
!     ---------------------------------------------------------------

      iShP2RS(:,1:NumAB)=0

      Do iSym = 1,nSym
         If (nAB(iSym) .gt. 0) Then
            iAB1 = iiBstRSh(iSym,iShlAB,iLoc) + 1
            iAB2 = iAB1 + nnBstRSh(iSym,iShlAB,iLoc) - 1
            Do iAB = iAB1,iAB2
               jAB = IndRed(iiBstR(iSym,iLoc)+iAB,iLoc) ! addr in 1st rs
               kAB = IndRed(jAB,1) ! addr in full shell pair
#if defined (_DEBUGPRINT_)
               nErr = 0
               If (IndRSh(jAB).ne.iSP2F(mySP(iShlAB))) Then
                  Write(Lupri,*) SecNam,': inconsistent shell pairs!'
                  Write(Lupri,*) SecNam,': from input: ',               &
     &                           iSP2F(mySP(iShlAB)),                   &
     &                           '  from IndRsh: ',IndRSh(jAB)
                  nErr = nErr + 1
               End If
               If (kAB.lt.1 .or. kAB.gt.NumAB) Then
                  Write(Lupri,*) SecNam,': shell pair address error!'
                  Write(Lupri,*) SecNam,': kAB = ',kAB
                  Write(Lupri,*) SecNam,': min and max allowed: 1 ',    &
     &                           NumAB
                  nErr = nErr + 1
               End If
               If (nErr .ne. 0) Then
                  Write(Lupri,*) SecNam,': Shell A, B, AB: ',           &
     &                           iShlA,iShlB,iShlAB
                  Write(Lupri,*) SecNam,': iLoc: ',iLoc
                  Write(Lupri,*) SecNam,': symmetry block: ',iSym
                  Write(Lupri,*) SecNam,': red. set address, ',         &
     &                           'first and current: ',                 &
     &                           jAB,iiBstR(iSym,iLoc)+iAB
                  Call Cho_Quit('Error detected in '//SecNam,104)
               End If
#endif
               iShP2RS(1,kAB)   = iAB
               iShP2RS(2,kAB) = iSym
            End Do
         End If
      End Do

!     Set return code 0: all ok!
!     --------------------------

      irc = 0

      End
