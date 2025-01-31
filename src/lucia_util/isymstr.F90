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

integer function ISYMSTR(ISYM,NSTR)
! Symmetry of product of NSTR string symmetries
!
! works currently only for D2H and subgroups
!
! Jeppe Olsen, 1998

use symmetry_info, only: MULTD2H => Mul

implicit none
! Input
integer ISYM(*), NSTR
integer IISYM, JSTR

if (NSTR == 0) then
  IISYM = 1
else
  IISYM = ISYM(1)
  do JSTR=2,NSTR
    IISYM = MULTD2H(IISYM,ISYM(JSTR))
  end do
end if

ISYMSTR = IISYM

end function ISYMSTR
