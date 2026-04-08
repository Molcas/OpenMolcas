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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ABSTR_TO_ORDSTR(IA_OC,IB_OC,NAEL,NBEL,IDET_OC,IDET_SP,SGN)
! An alpha string (IA) and a betastring (IB) is given.
! Combine these two strings to give an determinant with
! orbitals in ascending order. For doubly occupied orbitals
! the alphaorbital is given first.
! The output is given as IDET_OC : Orbital occupation (configuration)
!                        IDET_SP : Spin projections
!
! The phase required to change IA IB into IDET is computed as SGN
!
! Jeppe Olsen, November 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NAEL, IA_OC(NAEL), NBEL, IB_OC(NBEL)
integer(kind=iwp), intent(out) :: IDET_OC(NAEL+NBEL), IDET_SP(NAEL+NBEL), SGN
integer(kind=iwp) :: NEXT_AL, NEXT_BE, NEXT_EL

NEXT_AL = 1
NEXT_BE = 1
SGN = 1
! Loop over next electron in outputstring
do NEXT_EL=1,NAEL+NBEL
  if ((NEXT_AL <= NAEL) .and. (NEXT_BE <= NBEL)) then

    if (IA_OC(NEXT_AL) <= IB_OC(NEXT_BE)) then
      ! Next electron is alpha electron
      IDET_OC(NEXT_EL) = IA_OC(NEXT_AL)
      IDET_SP(NEXT_EL) = 1
      NEXT_AL = NEXT_AL+1
    else
      ! Next electron is beta electron
      IDET_OC(NEXT_EL) = IB_OC(NEXT_BE)
      IDET_SP(NEXT_EL) = -1
      NEXT_BE = NEXT_BE+1
      SGN = SGN*(-1)**(NAEL-NEXT_AL+1)
    end if
  else if (NEXT_BE > NBEL) then
    ! Next electron is alpha electron
    IDET_OC(NEXT_EL) = IA_OC(NEXT_AL)
    IDET_SP(NEXT_EL) = 1
    NEXT_AL = NEXT_AL+1
  else if (NEXT_AL > NAEL) then
    ! Next electron is beta electron
    IDET_OC(NEXT_EL) = IB_OC(NEXT_BE)
    IDET_SP(NEXT_EL) = -1
    NEXT_BE = NEXT_BE+1
    SGN = SGN*(-1)**(NAEL-NEXT_AL+1)
  end if
end do
! End of loop over electrons in outputlist

#ifdef _DEBUGPRINT_
write(u6,*) ' ABSTR to ORDSTR :'
write(u6,*) ' ================='
write(u6,*) ' Input alpha and beta strings'
call IWRTMA(IA_OC,1,NAEL,1,NAEL)
call IWRTMA(IB_OC,1,NBEL,1,NBEL)
write(u6,*) ' Configuration'
call IWRTMA(IDET_OC,1,NAEL+NBEL,1,NAEL+NBEL)
write(u6,*) ' Spin projections'
call IWRTMA(IDET_SP,1,NAEL+NBEL,1,NAEL+NBEL)
#endif

end subroutine ABSTR_TO_ORDSTR
