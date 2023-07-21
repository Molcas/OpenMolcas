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

subroutine CHO_OPEN(LUNIT,FNAME)
!
! Purpose: open sequential unformatted fortran file.

implicit none
integer LUNIT
character*(*) FNAME
integer LOCUNT, ISEED
integer ISFREEUNIT

if ((LUNIT < 1) .or. (LUNIT > 99)) then
  LOCUNT = 7
else
  LOCUNT = LUNIT
end if

ISEED = LOCUNT
LOCUNT = ISFREEUNIT(ISEED)
call molcas_binaryopen_vanilla(Locunt,Fname)
!open(LOCUNT,FILE=FNAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
LUNIT = LOCUNT

end subroutine CHO_OPEN
