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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

subroutine ZNELFSPGP(NTESTG)
! Generate for each supergroup the number of electrons in each active
! orbital space and store in NELFSPGP
!
! Jeppe Olsen, July 1995

use lucia_data, only: IBSPGPFTP, ISPGPFTP, MXPNGAS, NELFGP, NELFSPGP, NGAS, NSPGPFTP, NSTTP
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NTESTG
integer(kind=iwp) :: IBSPGP, IGAS, ISPGP, ITP, NSPGP, NTEST, NTESTL

NTESTL = 0
NTEST = max(NTESTG,NTESTL)

do ITP=1,NSTTP
  NSPGP = NSPGPFTP(ITP)
  IBSPGP = IBSPGPFTP(ITP)
  do ISPGP=IBSPGP,IBSPGP+NSPGP-1
    do IGAS=1,NGAS
      NELFSPGP(IGAS,ISPGP) = NELFGP(ISPGPFTP(IGAS,ISPGP))
    end do
  end do
end do

if (NTEST >= 10) then
  write(u6,*) ' Distribution of electrons in Active spaces'
  do ITP=1,NSTTP
    write(u6,*) ' String type ',ITP
    write(u6,*) ' Row : active space, Column: supergroup'
    NSPGP = NSPGPFTP(ITP)
    IBSPGP = IBSPGPFTP(ITP)
    call IWRTMA(NELFSPGP(1,IBSPGP),NGAS,NSPGP,MXPNGAS,NSPGP)
  end do
end if

end subroutine ZNELFSPGP
