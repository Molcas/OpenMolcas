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

subroutine GTSTTPS(IEL1,IEL3,NEL1,NEL3,NTYP,ITYP,IWAY)
! ITYP : type of strings with IEL1,IEL3 electrons
!
! IWAY = 1 : IEL1, IEL3 known, find ITYP
! IWAY = 2 : ITYP known, find IEL1, IEL3

implicit real*8(A-H,O-Z)
dimension NEL1(*), NEL3(*)

if (IWAY == 1) then
  ITYP = -1
  do IITYP=1,NTYP
    if ((IEL1 == NEL1(IITYP)) .and. (IEL3 == NEL3(IITYP))) ITYP = IITYP
  end do

  !if (ITYP == -1) then
  !  write(6,*) ' Error in GSTTPS'
  !  write(6,*) ' Error : Type could not be identified'
  !  write(6,*) ' Error : IEL1 IEL3 : ',IEL1,IEL3
  !  write(6,*) ' I am going to STOP'
  !  stop 'GSTTPS'
  !endif
else if (IWAY == 2) then
  IEL1 = NEL1(ITYP)
  IEL3 = NEL3(ITYP)
end if

NTEST = 0
if (NTEST >= 100) write(6,'(A,5I4)') ' GSTTPS : IWAY IEL1 IEL3 ITYP ',IWAY,IEL1,IEL3,ITYP

return

end subroutine GTSTTPS
