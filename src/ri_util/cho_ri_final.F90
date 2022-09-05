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
