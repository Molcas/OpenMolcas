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

subroutine MorsWrite(IMORS,STRING)

use Definitions, only: u6

implicit none
character(len=*) STRING
integer I, IB, IBIT, IMORS

if (IMORS < 0) goto 99
IB = IMORS
do I=1,len(STRING)
  STRING(I:I) = '0'
  IBIT = mod(IB,2)
  IB = IB/2
  if (IBIT == 1) STRING(I:I) = '1'
end do
if (IB > 0) then
  do I=1,len(STRING)
    STRING(I:I) = '*'
  end do
end if

return
99 continue
write(u6,*) ' MorsWrite: Bad IMORS=',IMORS
call ABEND()

end subroutine MorsWrite
