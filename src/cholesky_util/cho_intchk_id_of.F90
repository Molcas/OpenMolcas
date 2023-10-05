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

subroutine CHO_INTCHK_ID_OF(LABEL,ID,IOPT)
!
! Purpose: for minimal integral checking,
!          IOPT=-1 : return label corresponding to id ID.
!          else    : return index of shell quadruple corresponding
!                    to check label LABEL.

use Definitions, only: iwp

implicit none
character(len=8), intent(inout) :: LABEL
integer(kind=iwp), intent(inout) :: ID
integer(kind=iwp), intent(in) :: IOPT
integer(kind=iwp), parameter :: NTABLE = 12
character(len=*), parameter :: TABLE(NTABLE) = ['EXCL RS1','MAX|XRS1','MIN|XRS1','NEG DIAG','MAX|NEG ','MIN|NEG ','NEG->ZER', &
                                                'MAX|NEGZ','MIN|NEGZ','MAX DIAG','MIN DIAG','MAX|MIN ']
integer(kind=iwp), external :: CHO_TABIND

if (IOPT == -1) then
  if ((ID < 1) .or. (ID > NTABLE)) then
    LABEL = 'UNKNOWN '
  else
    LABEL = TABLE(ID)
  end if
else
  ID = CHO_TABIND(TABLE,8,NTABLE,' ',0,0,LABEL)
end if

end subroutine CHO_INTCHK_ID_OF
