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
      Subroutine bino(lmax)
      use Constants, only: Zero, One
      use welcom, only: binom
      Implicit None
      Integer lMax

      Integer i, j
!
      binom(:,:)=Zero
      binom(0,0)=One
      if(lmax.eq.0) Return
      Do i=1,lmax
         Do j=0,i
            binom(i,j)=binom(i-1,j-1)+binom(i-1,j)
         End Do
      End Do
!
      Return
      End Subroutine bino
