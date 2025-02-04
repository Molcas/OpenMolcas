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

subroutine ZSPGPIB(NSTSO,ISTSO,NSPGP,NSMST)
! Offset for supergroups of strings with given sym.
! Each supergroup is given start address 1
!
! Jeppe Olsen, Still summer of 95

use Definitions, only: u6

implicit real*8(A-H,O-Z)
! Input
integer NSTSO(NSMST,NSPGP)
! Output
integer ISTSO(NSMST,NSPGP)

do ISPGP=1,NSPGP
  ISTSO(1,ISPGP) = 1
  do ISMST=2,NSMST
    ISTSO(ISMST,ISPGP) = ISTSO(ISMST-1,ISPGP)+NSTSO(ISMST,ISPGP)
  end do
end do

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Output from ZSPGPIB'
  write(u6,*) ' ==================='
  write(u6,*)
  call IWRTMA(ISTSO,NSMST,NSPGP,NSMST,NSPGP)
end if

end subroutine ZSPGPIB
