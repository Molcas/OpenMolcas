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

module BasisMode

integer, parameter :: Valence_Mode = 0, Auxiliary_Mode = 1, Fragment_Mode = 2, With_Auxiliary_Mode = 3, With_Fragment_Mode = 4, &
                      All_Mode = 5
integer :: Basis_Mode, kCnttp, lCnttp
logical :: Atomic

end module BasisMode
