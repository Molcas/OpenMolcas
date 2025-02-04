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

subroutine MXMNOC_SPGP(MINEL,MAXEL,NORBTP,NORBFTP,NELFTP,NTESTG)
! Construct accumulated MAX and MIN arrays for a GAS supergroup

use Definitions, only: u6

implicit real*8(A-H,O-Z)
! Output
dimension MINEL(*), MAXEL(*)
! Input
integer NORBFTP(*), NELFTP(*)

! Some dummy initializations
IORB_START = 1 ! jwk-cleanup

NTESTL = 0
NTEST = max(NTESTG,NTESTL)

if (NTEST >= 100) then
  write(u6,*)
  write(u6,*) ' ==========='
  write(u6,*) ' MXMNOC_SPGP'
  write(u6,*) ' ==========='
  write(u6,*)
  !write(u6,*) ' NORBFTP :'
  !call IWRTMA(NORBFTP,1,NORBTP,1,NORBTP)
end if

do IORBTP=1,NORBTP
  ! Max and min at start of this type and at end of this type
  if (IORBTP == 1) then
    IORB_START = 1
    IORB_END = NORBFTP(1)
    NEL_START = 0
    NEL_END = NELFTP(1)
  else
    IORB_START = IORB_START+NORBFTP(IORBTP-1)
    IORB_END = IORB_START+NORBFTP(IORBTP)-1
    NEL_START = NEL_END
    NEL_END = NEL_START+NELFTP(IORBTP)
  end if
  if (NTEST >= 1000) then
    write(u6,*) ' IORBTP,IORB_START-IORB_END,NEL_START,NEL_END'
    write(u6,*) IORBTP,IORB_START-IORB_END,NEL_START,NEL_END
  end if

  do IORB=IORB_START,IORB_END
    MAXEL(IORB) = min(IORB,NEL_END)
    MINEL(IORB) = NEL_START
    if (NEL_END-MINEL(IORB) > IORB_END-IORB) MINEL(IORB) = NEL_END-(IORB_END-IORB)
  end do
end do

if (NTEST >= 100) then
  NORB = IELSUM(NORBFTP,NORBTP)
  write(u6,*) ' MINEL :'
  call IWRTMA(MINEL,1,NORB,1,NORB)
  write(u6,*) ' MAXEL :'
  call IWRTMA(MAXEL,1,NORB,1,NORB)
end if

end subroutine MXMNOC_SPGP
