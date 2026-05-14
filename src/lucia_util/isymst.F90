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

!#define _DEBUGPRINT_
function ISYMST(STRING,NEL)
! Symmmetry of string, D2H version

use Symmetry_Info, only: Mul
use lucia_data, only: ISMFTO
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: ISYMST
integer(kind=iwp), intent(in) :: NEL, STRING(NEL)
integer(kind=iwp) :: IEL

ISYMST = 1
do IEL=1,NEL
  ISYMST = Mul(ISYMST,ISMFTO(STRING(IEL)))
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' ISYMST, String and symmetry'
call IWRTMA(STRING,1,NEL,1,NEL)
write(u6,*) ISYMST
#endif

end function ISYMST
