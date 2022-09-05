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
      Logical Function Rsv_Tsk2(id,kls)
      use Tsk2
#include "stdalloc.fh"
      Logical, External :: Rsv_Tsk
!
      If (iOpt.eq.0) Then
         Rsv_Tsk2=Rsv_Tsk(id,kls)
      Else If (iOpt.eq.1) Then
         Rsv_Tsk2=.True.
         If (iRsv.gt.nTask) Then
            Rsv_Tsk2=.False.
         Else
            kls=TskList(iRsv)
            iRsv=iRsv+1
            If (kls.le.0) Rsv_Tsk2=.False.
            If (kls.gt.nTask) Rsv_Tsk2=.False.
         End If
      Else
         Rsv_Tsk2=.False.
         Call WarningMessage(2,'Error in Rsv_Tsk2')
         Write (6,*) 'Rsv_Tsk2: illegal iOpt value!'
         Call Abend()
      End If
!
      Return
      End
