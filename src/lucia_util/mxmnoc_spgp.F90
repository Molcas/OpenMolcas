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
subroutine MXMNOC_SPGP(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP)
! Construct accumulated MAX and MIN arrays for a GAS supergroup

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: MINEL(*), MAXEL(*)
integer(kind=iwp), intent(in) :: NORBTP, NORBFTP(NORBTP), NELFTP(NORBTP)
integer(kind=iwp) :: IORB, IORB_END, IORB_START, IORBTP, NEL_END, NEL_START
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NORB

write(u6,*)
write(u6,*) ' ==========='
write(u6,*) ' MXMNOC_SPGP'
write(u6,*) ' ==========='
write(u6,*)
!write(u6,*) ' NORBFTP :'
!call IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
#endif

IORB_START = 1
IORB_END = NORBFTP(1)
NEL_START = 0
NEL_END = NELFTP(1)
do IORBTP=1,NORBTP
  ! Max and min at start of this type and at end of this type
  if (IORBTP > 1) then
    IORB_START = IORB_START+NORBFTP(IORBTP-1)
    IORB_END = IORB_START+NORBFTP(IORBTP)-1
    NEL_START = NEL_END
    NEL_END = NEL_START+NELFTP(IORBTP)
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) ' IORBTP,IORB_START-IORB_END,NEL_START,NEL_END'
  write(u6,*) IORBTP,IORB_START-IORB_END,NEL_START,NEL_END
# endif

  do IORB=IORB_START,IORB_END
    MAXEL(IORB) = min(IORB,NEL_END)
    MINEL(IORB) = NEL_START
    if (NEL_END-MINEL(IORB) > IORB_END-IORB) MINEL(IORB) = NEL_END-(IORB_END-IORB)
  end do
end do

#ifdef _DEBUGPRINT_
NORB = sum(NORBFTP(1:NORBTP))
write(u6,*) ' MINEL :'
call IWRTMA(MINEL,1,NORB,1,NORB)
write(u6,*) ' MAXEL :'
call IWRTMA(MAXEL,1,NORB,1,NORB)
#endif

end subroutine MXMNOC_SPGP
