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
      Subroutine BasisConsistency(LuWr,iErr)
      Implicit Integer (i-n)
      Implicit Real*8 (a-h,o-z)
#include "g_zmatconv.fh"

      iErr = 0
      Do i = 1, 100
        If (BasReq(i) .AND. .NOT.BasAva(i)) GoTo 9900
      EndDo
      Return

9900  iErr = 1
      Write(LuWr,*) ' [BasisConsistency]: Atom NA=',i,' requires BS'
      Return
      End
