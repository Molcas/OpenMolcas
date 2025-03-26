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

subroutine RSMXMN_MCLR(MAXEL,MINEL,NORB1,NORB2,NORB3,NEL,MIN1,MAX1,MIN3,MAX3)
! Construct accumulated MAX and MIN arrays for a RAS set of strings

implicit real*8(A-H,O-Z)
dimension MINEL(*), MAXEL(*)

! accumulated max and min in each of the three spaces
! (required max and min at final orbital in each space)
MIN1A = MIN1
MAX1A = MAX1

MIN2A = NEL-MAX3
MAX2A = NEL-MIN3

MIN3A = NEL
MAX3A = NEL

NORB = NORB1+NORB2+NORB3
do IORB=1,NORB
  if (IORB <= NORB1) then
    MINEL(IORB) = max(MIN1A+IORB-NORB1,0)
    MAXEL(IORB) = min(IORB,MAX1A)
  else if ((NORB1 < IORB) .and. (IORB <= (NORB1+NORB2))) then
    MINEL(IORB) = max(MIN2A+IORB-NORB1-NORB2,0)
    if (NORB1 > 0) MINEL(IORB) = max(MINEL(IORB),MINEL(NORB1))
    MAXEL(IORB) = min(IORB,MAX2A)
  else if (IORB > NORB1+NORB2) then
    MINEL(IORB) = max(MIN3A+IORB-NORB,0)
    if (NORB1+NORB2 > 0) MINEL(IORB) = max(MINEL(IORB),MINEL(NORB1+NORB2))
    MAXEL(IORB) = min(IORB,MAX3A)
  end if
end do

return

end subroutine RSMXMN_MCLR
