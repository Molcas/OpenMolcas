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
      Subroutine Off_Diagonal(B1,nB,iB1s,iB1e,B2,iB2s,iB2e)
      Implicit Real*8 (a-h,o-z)
      Real*8 B1(nB,iB1s:iB1e), B2(nB,iB2s:iB2e)
!
      Do j = iB2s, iB2e
         Do i = iB1s, iB1e
            B1(j,i)=B2(i,j)
         End Do
      End Do
!
      Return
      End
