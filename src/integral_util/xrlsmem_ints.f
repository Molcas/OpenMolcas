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
      Subroutine xRlsMem_Ints
      implicit real*8 (a-h,o-z)
#include "status.fh"
!
      If (XMem_Status.eq.InActive) Then
!        Write (6,*)
!        Write (6,*) 'xRlsMem_Ints:',
!    &               'No external scratch handling to deactivate!'
!        Write (6,*)
      Else
!        Write (6,*)
!        Write (6,*) 'xRlsMem_Ints:','External allocation deactivate!'
!        Write (6,*)
         XMem_Status=InActive
         Call RlsMem_Ints()
!        Write (6,*) 'xRlsMem_Ints:','External allocation released!'
      End If
!
      Return
      End
