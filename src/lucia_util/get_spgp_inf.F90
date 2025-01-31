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

use lucia_data, only: NGAS
use lucia_data, only: IBSPGPFTP, ISPGPFTP

implicit none
integer ISPGP, ITP
integer IGRP(*)
integer NTEST, ISPGPABS

NTEST = 0
! Absolute group number
!write(6,*) ' GET_SPGP_INF : ISPGP, ITP',ISPGP,ITP
ISPGPABS = ISPGP+IBSPGPFTP(ITP)-1
call ICOPVE(ISPGPFTP(1,ISPGPABS),IGRP,NGAS)

if (NTEST >= 100) then
  write(6,*) ' GET_SPGP_INF : ISPGP ITP ISPGPABS',ISPGP,ITP,ISPGPABS
  write(6,*) ' Groups of supergroups'
  call IWRTMA(IGRP,1,NGAS,1,NGAS)
end if

end subroutine GET_SPGP_INF
