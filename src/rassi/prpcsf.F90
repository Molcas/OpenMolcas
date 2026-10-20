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

subroutine PRPCSF(IOPEN,NCPL,ICOUP)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOPEN, NCPL, ICOUP(IOPEN,NCPL)
integer(kind=iwp) :: I, J, N
character(len=24) :: FRMT
character, parameter :: CPLSMB(0:1) = ['d','u']

if ((IOPEN < 0) .or. (NCPL < 0)) then
  call WarningMessage(2,'Program bug: Erroneous call to PRPCSF.')
  write(u6,*) 'PRPCSF error: Wrong arguments.'
  write(u6,*) 'PRPCSF: IOPEN=',IOPEN
  write(u6,*) 'PRPCSF: NCPL =',NCPL
  call ABEND()
end if
if ((IOPEN == 0) .or. (NCPL == 0)) then
  call WarningMessage(1,'Program bug? Strange call to PRPCSF.')
  write(u6,*) 'PRPCSF warning: Strange arguments.'
  write(u6,*) 'PRPCSF: IOPEN=',IOPEN
  write(u6,*) 'PRPCSF: NCPL =',NCPL
else
  N = 80/(7+IOPEN)
  write(FRMT,'(A1,I2,A,I2,A4)') '(',N,'(1X,I5,1X,',IOPEN,'A1))'
  write(u6,FRMT) (I,(CPLSMB(ICOUP(J,I)),J=1,IOPEN),I=1,NCPL)
end if

end subroutine PRPCSF
