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

!#define _DEBUGPRINT_
subroutine GTSTTPS_1(IEL1,IEL3,NEL1,NEL3,NTYP,ITYP)
! ITYP : type of strings with IEL1,IEL3 electrons
!
! IEL1, IEL3 known, find ITYP

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IEL1, IEL3, NEL1(*), NEL3(*), NTYP
integer(kind=iwp), intent(out) :: ITYP
integer(kind=iwp) :: IITYP

ITYP = -1
do IITYP=1,NTYP
  if ((IEL1 == NEL1(IITYP)) .and. (IEL3 == NEL3(IITYP))) ITYP = IITYP
end do

!if (ITYP == -1) then
!  write(u6,*) ' Error in GTSTTPS_1'
!  write(u6,*) ' Error : Type could not be identified'
!  write(u6,*) ' Error : IEL1 IEL3 : ',IEL1,IEL3
!  write(u6,*) ' I am going to STOP'
!  stop 'GSTTPS'
!endif

#ifdef _DEBUGPRINT_
write(u6,'(A,5I4)') ' GTSTTPS_1 : IEL1 IEL3 ITYP ',IEL1,IEL3,ITYP
#endif

end subroutine GTSTTPS_1
