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
!     This defines the highest angular momentum quantum number which
!     Seward will be able to treat.
Module Define_af
#ifndef _DEMO_
      Integer, parameter:: iTabMx=15
      Character(LEN=1), Parameter:: AngTp(0:iTabMx) = ['s','p','d','f','g',                          &
     &                                                 'h','i','k','l','m',                          &
     &                                                 'n','o','q','r','t',                          &
     &                                                 'u']
#else
      Integer, parameter:: iTabMx=3
      Character(LEN=1), Parameter:: AngTp(0:iTabMx) = ['s','p','d']
#endif
End Module Define_af