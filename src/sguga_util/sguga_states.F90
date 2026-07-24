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

module SGUGA_States

use sguga, only: SGStruct, EXStruct, CIStruct
use definitions, only: iwp

Type (SGStruct) :: SGS(3)
Type (CIStruct) :: CIS(3)
Type (EXStruct) :: EXS(3)
Logical(kind=iwp) :: State_is_used(2)=[.False.,.False.]

Public:: SGS, CIS, EXS, State_is_used

end module SGUGA_States
