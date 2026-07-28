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

use sguga, only: CIStruct, EXStruct, SGStruct
use Definitions, only: iwp

type(SGStruct) :: SGS(3)
type(CIStruct) :: CIS(3)
type(EXStruct) :: EXS(3)

public :: CIS, EXS, SGS

end module SGUGA_States
