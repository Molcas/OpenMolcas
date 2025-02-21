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
! Copyright (C) 1990, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
function ISTRNM(IOCC,NORB,NEL,Z,NEWORD,IREORD)
! Address of string IOCC
!
! version of Winter 1990, Jeppe Olsen

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: ISTRNM
integer(kind=iwp), intent(in) :: NEL, IOCC(NEL), NORB, Z(NORB,NEL), NEWORD(*), IREORD
integer(kind=iwp) :: I, IZ

IZ = 1
do I=1,NEL
  IZ = IZ+Z(IOCC(I),I)
end do

if (IREORD == 0) then
  ISTRNM = IZ
else
  ISTRNM = NEWORD(IZ)
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' STRING'
call IWRTMA(IOCC,1,NEL,1,NEL)
write(u6,*) ' Z matrix'
call IWRTMA(Z,NORB,NEL,NORB,NEL)
!write(u6,*) ' First two elements of reorder array'
!call IWRTMA(NEWORD,1,2,1,2)
write(u6,*) ' ADDRESS OF STRING ',ISTRNM
write(u6,*) ' REV LEX number : ',IZ
#endif

end function ISTRNM
