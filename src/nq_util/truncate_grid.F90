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
      Subroutine Truncate_Grid(R,nR,nR_Eff,Radius_Max)
      Implicit Real*8 (a-h,o-z)
      Real*8 R(2,nR)
!
      nTmp=nR_Eff
      Do i = 1, nTmp
         If (R(1,i).gt.Radius_Max) Then
             nR_Eff=i-1
             Go To 99
         End If
      End Do
 99   Continue
!
      Return
      End
