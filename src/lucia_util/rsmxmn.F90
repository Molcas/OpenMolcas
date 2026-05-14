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
subroutine RSMXMN(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,MN1,MX1,MN3,MX3)
! Construct accumulated MAX and MIN arrays for a RAS set of strings

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NORB1, NORB2, NORB3, NEL, MN1, MX1, MN3, MX3
integer(kind=iwp), intent(out) :: MAXEL(NORB1+NORB2+NORB3), MINEL(NORB1+NORB2+NORB3)
integer(kind=iwp) :: IORB, MAX1A, MAX2A, MAX3A, MIN1A, MIN2A, MIN3A, NORB

NORB = NORB1+NORB2+NORB3
! accumulated max and min in each of the three spaces
! (required max and min at final orbital in each space)
!OLD MIN1A = MN1
MIN1A = max(MN1,NEL-MX3-NORB2)
MAX1A = MX1

MIN2A = NEL-MX3
MAX2A = NEL-MN3

MIN3A = NEL
MAX3A = NEL

do IORB=1,NORB
  if (IORB <= NORB1) then
    MINEL(IORB) = max(MIN1A+IORB-NORB1,0)
    MAXEL(IORB) = min(IORB,MAX1A)
  else if (IORB <= (NORB1+NORB2)) then
    MINEL(IORB) = max(MIN2A+IORB-NORB1-NORB2,0)
    if (NORB1 > 0) MINEL(IORB) = max(MINEL(IORB),MINEL(NORB1))
    MAXEL(IORB) = min(IORB,MAX2A)
  else
    MINEL(IORB) = max(MIN3A+IORB-NORB,0)
    if (NORB1+NORB2 > 0) MINEL(IORB) = max(MINEL(IORB),MINEL(NORB1+NORB2))
    MAXEL(IORB) = min(IORB,MAX3A)
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from RSMXMN'
write(u6,*) ' =================='
write(u6,*) ' MINEL :'
call IWRTMA(MINEL,1,NORB,1,NORB)
write(u6,*) ' MAXEL :'
call IWRTMA(MAXEL,1,NORB,1,NORB)
#endif

end subroutine RSMXMN
