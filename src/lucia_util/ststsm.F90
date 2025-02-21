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
subroutine STSTSM(STSTSX,STSTDX,NSMST)
! construct STSTSX and STSTDX giving
! symmetry of sx (dx) connecting two given string symmetries

use Symmetry_Info, only: Mul
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSMST
integer(kind=iwp), intent(out) :: STSTSX(NSMST,NSMST), STSTDX(NSMST,NSMST)
integer(kind=iwp) :: ILSTSM, IRSTSM

do ILSTSM=1,NSMST
  do IRSTSM=1,NSMST
    STSTSX(ILSTSM,IRSTSM) = Mul(IRSTSM,ILSTSM)
    STSTDX(ILSTSM,IRSTSM) = Mul(IRSTSM,ILSTSM)
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' STSTSM : STSTSX, STSTDX'
call IWRTMA(STSTSX,NSMST,NSMST,NSMST,NSMST)
call IWRTMA(STSTDX,NSMST,NSMST,NSMST,NSMST)
#endif

end subroutine STSTSM
