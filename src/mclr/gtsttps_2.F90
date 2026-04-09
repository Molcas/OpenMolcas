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
subroutine GTSTTPS_2(IEL1,IEL3,NEL1,NEL3,ITYP)
! ITYP : type of strings with IEL1,IEL3 electrons
!
! ITYP known, find IEL1, IEL3

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(out) :: IEL1, IEL3
integer(kind=iwp), intent(in) :: NEL1(*), NEL3(*), ITYP

IEL1 = NEL1(ITYP)
IEL3 = NEL3(ITYP)

#ifdef _DEBUGPRINT_
write(u6,'(A,5I4)') ' GTSTTPS : IEL1 IEL3 ITYP ',IEL1,IEL3,ITYP
#endif

end subroutine GTSTTPS_2
