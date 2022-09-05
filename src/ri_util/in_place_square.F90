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
      Subroutine In_place_Square(Buff,nBuff)
      Implicit Real*8 (a-h,o-z)
      Real*8 Buff(nBuff,nBuff)
!
!     Call RecPrt('Buff',' ',Buff,nBuff,nBuff)
      Do j = 1, nBuff
         Do i = 1, j-1
            Buff(j,i)=Buff(i,j)
         End Do
      End Do
!     Call RecPrt('Buff',' ',Buff,nBuff,nBuff)
!     Write (6,'(10F10.3)') (Buff(i,i),i=1,nBuff)
!
      Return
      End
