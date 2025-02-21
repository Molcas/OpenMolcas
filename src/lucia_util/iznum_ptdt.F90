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

!#define _DEBUGPRINT_
function IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
! Address of prototype determinant IAB
! alpha occupation is used to define lex address
!
! Jeppe Olsen, Dec. 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IZNUM_PTDT
integer(kind=iwp), intent(in) :: NOPEN, IAB(NOPEN), NALPHA, Z(NOPEN,NALPHA), NEWORD(*), IREORD
integer(kind=iwp) :: I, IALPHA, IZ

#ifdef _DEBUGPRINT_
write(u6,*) ' IZNUM_PTDT, NOPEN, NALPHA = ',NOPEN,NALPHA
write(u6,*) ' Input Z- matrix'
call IWRTMA(Z,NOPEN,NALPHA,NOPEN,NALPHA)
#endif

IZ = 1
IALPHA = 0
do I=1,NOPEN
  if (IAB(I) > 0) then
    IALPHA = IALPHA+1
    IZ = IZ+Z(I,IALPHA)
  end if
end do

!write(u6,*) ' IZ = ',IZ
if (IREORD == 0) then
  IZNUM_PTDT = IZ
else
  IZNUM_PTDT = NEWORD(IZ)
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from IZNUM_PTDT'
write(u6,*) ' Prototype determinant'
call IWRTMA(IAB,1,NOPEN,1,NOPEN)
write(u6,*) ' Address = ',IZNUM_PTDT
#endif

end function IZNUM_PTDT
