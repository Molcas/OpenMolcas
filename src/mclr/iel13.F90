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
! Copyright (C) 1991,1994, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine IEL13(MNRS1,MXRS1,MNRS3,MXRS3,NELEC,NOCTYP,IEL1,IEL3,IEL123)
! A type of strings contains NEL electrons and
!
! The number of electrons in RAS 1 is between MNRS1 AND MXRS1
! The number of electrons in RAS 3 is between MNRS3 AND MXRS3
!
! Find the number of electrons in RAS1 and RAS3 in each subtype
! and determine if subtype  is active
!
! Jeppe Olsen, Spring 1991
!              IEL123 added Jan 1994

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MNRS1, MXRS1, MNRS3, MXRS3, NELEC, NOCTYP
integer(kind=iwp), intent(out) :: IEL1(NOCTYP), IEL3(NOCTYP), IEL123(3,NOCTYP)
integer(kind=iwp) :: ITYP, KEL1, KEL3, NTYP

IEL1(:) = 0
IEL3(:) = 0
do KEL3=MNRS3,MXRS3
  do KEL1=MNRS1,MXRS1
    ITYP = (MXRS1-KEL1)*(MXRS3-MNRS3+1)+KEL3-MNRS3+1
    IEL1(ITYP) = KEL1
    IEL3(ITYP) = KEL3
  end do
end do
NTYP = (MXRS1-MNRS1+1)*(MXRS3-MNRS3+1)
IEL123(1,1:NTYP) = IEL1(1:NTYP)
IEL123(2,1:NTYP) = NELEC-IEL1(1:NTYP)-IEL3(1:NTYP)
IEL123(3,1:NTYP) = IEL3(1:NTYP)

#ifdef _DEBUGPRINT_
write(u6,*) ' =============='
write(u6,*) ' IEL13 speaking'
write(u6,*) ' =============='
write(u6,'(A,4I3)') ' IEL1 IEL3 for MNRS1 MXRS1 MNRS3 MXRS3 ',MNRS1,MXRS1,MNRS3,MXRS3
do ITYP=1,NOCTYP
  write(u6,'(3I3)') IEL1(ITYP),IEL3(ITYP)
end do

write(u6,*) ' IEL123 matrix'
call IWRTMA(IEL123,3,NOCTYP,3,NOCTYP)
#endif

return

end subroutine IEL13
