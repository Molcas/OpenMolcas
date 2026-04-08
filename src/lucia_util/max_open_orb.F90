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
subroutine MAX_OPEN_ORB(MAXOP,IOCLS,NGAS,NOCLS,NOBPT)
! Max number of open orbitals in occupation classes
!
! Jeppe Olsen, November 2001

use csfbas, only: maxop_lucia
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(out) :: MAXOP
integer(kind=iwp), intent(in) :: NGAS, NOCLS, IOCLS(NGAS,NOCLS), NOBPT(NGAS)
integer(kind=iwp) :: IGAS, JOCLS, MAXOP_IGAS, MAXOP_J, NEL, NORB

MAXOP = 0
!write(u6,*) ' NOCLS, NGAS = ',NOCLS,NGAS
do JOCLS=1,NOCLS
  MAXOP_J = 0
  do IGAS=1,NGAS
    NEL = IOCLS(IGAS,JOCLS)
    NORB = NOBPT(IGAS)
    !write(u6,*) ' IGAS, NEL, NORB = ',IGAS,NEL,NORB
    MAXOP_IGAS = min(NEL,2*NORB-NEL)
    MAXOP_J = MAXOP_J+MAXOP_IGAS
  end do
  MAXOP = max(MAXOP,MAXOP_J)
end do
maxop_lucia = maxop

#ifdef _DEBUGPRINT_
write(u6,*) ' Max number of unpaired orbitals = ',MAXOP
#endif

end subroutine MAX_OPEN_ORB
