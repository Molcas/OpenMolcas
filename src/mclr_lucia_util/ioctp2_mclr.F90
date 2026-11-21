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

function IOCTP2_MCLR(STRING,NEL,ITYP)
! Obtain occupation type for STRING.
! For forbidden strings a zero is returned
!
! New version allowing general set of strings

use Str_Info, only: MNRS1, MNRS3, MXRS1, MXRS3
use MCLR_Lucia_Data, only: NORB1, NORB2
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: IOCTP2_MCLR
integer(kind=iwp), intent(in) :: NEL, STRING(NEL), ITYP
integer(kind=iwp) :: iEl, iEL1, iEL3, iTyp2

if (ITYP <= 0) then
  write(u6,*) 'IOCTP2: ITYP <= 0'
  write(u6,*) 'ITYP=',ITYP
  call Abend()
end if
! Number of electrons in RAS1 and RAS 3
IEL1 = 0
IEL3 = 0
do IEL=1,NEL
  if (STRING(IEL) <= NORB1) IEL1 = IEL1+1
  if (NORB1+NORB2+1 <= STRING(IEL)) IEL3 = IEL3+1
end do
! Type
if (((IEL1 >= MNRS1(ITYP)) .and. (IEL1 <= MXRS1(ITYP))) .and. ((IEL3 >= MNRS3(ITYP)) .and. (IEL3 <= MXRS3(ITYP)))) then
  ITYP2 = (MXRS1(ITYP)-IEL1)*(MXRS3(ITYP)-MNRS3(ITYP)+1)+IEL3-MNRS3(ITYP)+1
else
  ITYP2 = 0
end if
IOCTP2_MCLR = ITYP2

end function IOCTP2_MCLR
