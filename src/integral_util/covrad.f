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
      Real*8 Function CovRad(i)
      use CovRad_Data, only: CovRad_
      Implicit None
      Integer i
!
      If (i.gt.86) Then
!        Write (Warning,'(2A)') 'CovRad: Warning i.gt.86!,;'//
!    &               'Guestimate of 2.70 au is used!'
!        Call WarningMessage(1,Warning)
         CovRad=2.70D0
      Else
         CovRad=CovRad_(i)
      End If
!
      Return
      End Function CovRad
