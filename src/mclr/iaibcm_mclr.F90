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
subroutine IAIBCM_MCLR(MNRS1,MXRS3,NOCTPA,NOCTPB,NEL1A,NEL3A,NEL1B,NEL3B,IOCOC)
! RAS allowed combinations of alpha and beta types
!
! =====
! Input
! =====
!
! NOCTPA : Number of alpha types
! NEL1A  : Number of electrons in RAS 1 for each alpha type
! NEL3A  : Number of electrons in RAS 3 for each alpha type
!
! NOCTPB : Number of beta types
! NEL1B  : Number of electrons in RAS 1 for each beta type
! NEL3B  : Number of electrons in RAS 3 for each beta type
!
! ======
! Output
! ======
! IOCOC(IATP,IBTP)  = 1 =>      allowed combination
! IOCOC(IATP,IBTP)  = 0 => not allowed combination

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MNRS1, MXRS3, NOCTPA, NOCTPB, NEL1A(NOCTPA), NEL3A(NOCTPA), NEL1B(NOCTPB), NEL3B(NOCTPB)
integer(kind=iwp), intent(out) :: IOCOC(NOCTPA,NOCTPB)
integer(kind=iwp) :: IAEL1, IAEL3, IATP, IBEL1, IBEL3, IBTP

IOCOC(:,:) = 0
do IATP=1,NOCTPA
  IAEL1 = NEL1A(IATP)
  IAEL3 = NEL3A(IATP)
  do IBTP=1,NOCTPB
    IBEL1 = NEL1B(IBTP)
    IBEL3 = NEL3B(IBTP)
    if ((IAEL1+IBEL1 >= MNRS1) .and. (IAEL3+IBEL3 <= MXRS3)) IOCOC(IATP,IBTP) = 1
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' Matrix giving allowed combinations of types'
write(u6,*)
call IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
#endif

return

end subroutine IAIBCM_MCLR
