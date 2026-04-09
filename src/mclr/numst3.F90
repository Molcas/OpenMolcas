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
function NUMST3(NEL,NORB1,NEL1MN,NEL1MX,NORB2,NORB3,NEL3MN,NEL3MX)
! Number of strings with NEL electrons that fulfill
!
! Between NEL1MN AND NEL1MX electrons in the first NORB1 orbitals
! Between NEL3MN AND NEL3MX electrons in the last  NORB3 orbitals

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: NUMST3
integer(kind=iwp), intent(in) :: NEL, NORB1, NEL1MN, NEL1MX, NORB2, NORB3, NEL3MN, NEL3MX
integer(kind=iwp) :: IEL1, IEL2, IEL3, IEL3MN, IEL3MX, NSTIN1, NSTINT, NSTRIN
integer(kind=iwp), external :: IBINOM

NSTRIN = 0

do IEL1=NEL1MN,min(NEL1MX,NORB1,NEL)
  NSTIN1 = IBINOM(NORB1,IEL1)
  IEL3MN = max(NEL3MN,NEL-(IEL1+NORB2))
  IEL3MX = min(NEL3MX,NEL-IEL1)
  do IEL3=IEL3MN,IEL3MX
    IEL2 = NEL-IEL1-IEL3
    NSTINT = NSTIN1*IBINOM(NORB2,IEL2)*IBINOM(NORB3,IEL3)
    NSTRIN = NSTRIN+NSTINT
  end do
end do
NUMST3 = NSTRIN

#ifdef _DEBUGPRINT_
write(u6,'(/A,I6)') '  Number of strings generated ... ',NSTRIN
#endif

return

end function NUMST3
