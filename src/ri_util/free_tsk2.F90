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
      Subroutine Free_Tsk2(id)
      use Tsk2
#include "stdalloc.fh"
!
      If (iOpt.eq.0) Then
         Call Free_Tsk(id)
      Else If (iOpt.eq.1) Then
         Call mma_deallocate(TskList)
         nTask=0
      Else
         Call WarningMessage(2,'Error in Free_Tsk2')
         Write (6,*) 'Free_Tsk2: illegal iOpt value!'
         Call Abend()
      End If
      iOpt=-1
!
      Return
      End
