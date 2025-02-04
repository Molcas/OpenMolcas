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

integer function ISYMS1(STRING,NEL)
! Symmmetry of string, D2H version

use symmetry_info, only: SYMPRO => Mul
use lucia_data, only: ISMFTO
use Definitions, only: u6

implicit none
! Specific input
integer NEL
integer STRING(*)
integer ISYM, IEL, NTEST

ISYM = 1
do IEL=1,NEL
  ISYM = SYMPRO(ISYM,ISMFTO(STRING(IEL)))
end do

ISYMS1 = ISYM

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' ISYMS1, String and symmetry'
  call IWRTMA(STRING,1,NEL,1,NEL)
  write(u6,*) ISYM
end if

end function ISYMS1
