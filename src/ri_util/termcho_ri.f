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
      SubRoutine TermCho_RI(irc,nVec_RI,l_nVec_RI)
      Implicit None
      Integer irc, l_nVec_RI
      Integer nVec_RI(l_nVec_RI) ! #RI vectors per irrep on this node

      irc = 0

C     Save number of vectors and other info on runfile.
C     -------------------------------------------------

      Call Cho_Final(.False.)
      Call Cho_RI_Final(irc,nVec_RI,l_nVec_RI)
      If (irc .ne. 0) Return

C     Close storage files.
C     --------------------

      Call Cho_P_OpenVR(2)

C     Deallocate index arrays.
C     ------------------------

      Call Cho_X_Dealloc(irc)
      If (irc .ne. 0) Return

C     More deallocations.
C     -------------------

      Call Cho_RI_XFree()

      End
      SubRoutine Cho_RI_XFree()
      Implicit None
#include "choptr2.fh"

      If (l_mySP .gt. 0) Then
         Call GetMem('mySP','Free','Inte',ip_mySP,l_mySP)
      End If

      End
      SubRoutine Cho_RI_Final(irc,nVec_RI,l_nVec_RI)
      Implicit None
      Integer irc, l_nVec_RI
      Integer nVec_RI(l_nVec_RI)
#include "cholesky.fh"

      If (l_nVec_RI .lt. nSym) Then
         irc = 1
         Return
      Else
         irc = 0
         Call Put_iArray('nVec_RI',nVec_RI,nSym)
      End If

      End
