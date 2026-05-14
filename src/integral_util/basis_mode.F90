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

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: Valence_Mode = 0, Auxiliary_Mode = 1, Fragment_Mode = 2, With_Auxiliary_Mode = 3, &
                                With_Fragment_Mode = 4, All_Mode = 5
integer(kind=iwp) :: Basis_Mode, kCnttp, lCnttp
logical(kind=iwp) :: Atomic

public :: All_Mode, Atomic, Auxiliary_Mode, Basis_Mode, Fragment_Mode, kCnttp, lCnttp, Valence_Mode, With_Auxiliary_Mode, &
          With_Fragment_Mode

end module BasisMode
