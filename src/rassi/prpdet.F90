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

subroutine PRPDET(IOPEN,ND,IDET)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOPEN, ND, IDET(IOPEN,ND)
integer(kind=iwp) :: I, J, N
character(len=24) :: FRMT
character, parameter :: SPNSMB(0:1) = ['b','a']

if ((IOPEN < 0) .or. (ND < 0)) then
  call WarningMessage(2,'Program bug: Erroneous call to PRPDET.')
  write(u6,*) 'PRPDET error: Wrong arguments.'
  write(u6,*) 'PRPDET: IOPEN=',IOPEN
  write(u6,*) 'PRPDET: ND =',ND
  call ABEND()
end if
if ((IOPEN == 0) .or. (ND == 0)) then
  call WarningMessage(1,'Program bug? Strange call to PRPDET.')
  write(u6,*) 'PRPDET warning: Strange arguments.'
  write(u6,*) 'PRPDET: IOPEN=',IOPEN
  write(u6,*) 'PRPDET: ND =',ND
else
  N = 80/(7+IOPEN)
  write(FRMT,'(A1,I2,A,I2,A4)') '(',N,'(1X,I5,1X,',IOPEN,'A1))'
  write(u6,FRMT) (I,(SPNSMB(IDET(J,I)),J=1,IOPEN),I=1,ND)
end if

end subroutine PRPDET
