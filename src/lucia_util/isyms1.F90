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

function ISYMS1(STRING,NEL)
! Symmmetry of string, D2H version

use Symmetry_Info, only: Mul
use lucia_data, only: ISMFTO
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ISYMS1
integer(kind=iwp) :: STRING(*), NEL
integer(kind=iwp) :: IEL, ISYM, NTEST

ISYM = 1
do IEL=1,NEL
  ISYM = Mul(ISYM,ISMFTO(STRING(IEL)))
end do

ISYMS1 = ISYM

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' ISYMS1, String and symmetry'
  call IWRTMA(STRING,1,NEL,1,NEL)
  write(u6,*) ISYM
end if

end function ISYMS1
