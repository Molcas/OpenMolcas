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
! Copyright (C) 1998, Jeppe Olsen                                      *
!***********************************************************************

function ISYMSTR(ISYM,NSTR)
! Symmetry of product of NSTR string symmetries
!
! works currently only for D2H and subgroups
!
! Jeppe Olsen, 1998

use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ISYMSTR
integer(kind=iwp), intent(in) :: NSTR, ISYM(NSTR)
integer(kind=iwp) :: JSTR

ISYMSTR = 1
do JSTR=1,NSTR
  ISYMSTR = Mul(ISYMSTR,ISYM(JSTR))
end do

end function ISYMSTR
