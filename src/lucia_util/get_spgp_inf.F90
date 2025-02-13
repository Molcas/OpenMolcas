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
! Copyright (C) 1997, Jeppe Olsen                                      *
!***********************************************************************

subroutine GET_SPGP_INF(ISPGP,ITP,IGRP)
! Obtain groups defining supergroup ISPGP,ITP
!
! Jeppe Olsen, November 97

use lucia_data, only: IBSPGPFTP, ISPGPFTP, NGAS
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ISPGP, ITP
integer(kind=iwp), intent(out) :: IGRP(NGAS)
integer(kind=iwp) :: ISPGPABS, NTEST

NTEST = 0
! Absolute group number
!write(u6,*) ' GET_SPGP_INF : ISPGP, ITP',ISPGP,ITP
ISPGPABS = ISPGP+IBSPGPFTP(ITP)-1
IGRP(:) = ISPGPFTP(1:NGAS,ISPGPABS)

if (NTEST >= 100) then
  write(u6,*) ' GET_SPGP_INF : ISPGP ITP ISPGPABS',ISPGP,ITP,ISPGPABS
  write(u6,*) ' Groups of supergroups'
  call IWRTMA(IGRP,1,NGAS,1,NGAS)
end if

end subroutine GET_SPGP_INF
