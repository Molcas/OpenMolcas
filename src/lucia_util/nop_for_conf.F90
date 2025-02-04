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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

function NOP_FOR_CONF(ICONF,NEL)
! A configuration is given as a nonstrict ascending sequence of occupied
! occupied orbitals. Find number of double occupied orbitals
!
! Jeppe Olsen, Nov. 2001

use Definitions, only: u6

implicit real*8(A-H,O-Z)
integer ICONF(NEL)

! Loop over electrons
NOPEN = 0
IEL = 1
1000 continue
if (IEL < NEL) then
  if (ICONF(IEL) /= ICONF(IEL+1)) then
    NOPEN = NOPEN+1
    IEL = IEL+1
  else if (ICONF(IEL) == ICONF(IEL+1)) then
    IEL = IEL+2
  end if
end if

if (IEL == NEL) then
  ! The last orbital is not identical to any later orbitals so
  NOPEN = NOPEN+1
  IEL = IEL+1
end if
if (IEL < NEL) goto 1000

NOP_FOR_CONF = NOPEN

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Configuration'
  call IWRTMA(ICONF,1,NEL,1,NEL)
  write(u6,*) ' Number of open orbitals = ',NOP_FOR_CONF
end if

end function NOP_FOR_CONF
