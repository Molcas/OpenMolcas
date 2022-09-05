!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2008, Roland Lindh                                     *
!***********************************************************************
      Subroutine ReNorm2(iCnttp)
      use Wrj12, only: iOffA
!
      Call ICopy(4*8,[0],0,iOffA,1)
      Do ire_do = 1, 2
!
         Call ReNorm2_(iCnttp)
!
      End Do
!
      Return
      End
