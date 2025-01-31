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

function IZNUM_PTDT(IAB,NOPEN,NALPHA,Z,NEWORD,IREORD)
! Address of prototype determinant IAB
! alpha occupation is used to define lex address
!
! Jeppe Olsen, Dec. 2001

implicit real*8(A-H,O-Z)
integer Z(NOPEN,NALPHA)
dimension IAB(*), NEWORD(*)

NTEST = 0
if (NTEST >= 100) then
  write(6,*) ' IZNUM_PTDT, NOPEN, NALPHA = ',NOPEN,NALPHA
  write(6,*) ' Input Z- matrix'
  call IWRTMA(Z,NOPEN,NALPHA,NOPEN,NALPHA)
end if

IZ = 1
IALPHA = 0
do I=1,NOPEN
  if (IAB(I) > 0) then
    IALPHA = IALPHA+1
    IZ = IZ+Z(I,IALPHA)
  end if
end do

!write(6,*) ' IZ = ',IZ
if (IREORD == 0) then
  IZNUM_PTDT = IZ
else
  IZNUM_PTDT = NEWORD(IZ)
end if

if (NTEST >= 100) then
  write(6,*) ' Output from IZNUM_PTDT'
  write(6,*) ' Prototype determinant'
  call IWRTMA(IAB,1,NOPEN,1,NOPEN)
  write(6,*) ' Address = ',IZNUM_PTDT
end if

end function IZNUM_PTDT
