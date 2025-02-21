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
function ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
! A configuration ICONF of NOCC_ORB orbitals are given
! ICONF(I) = IORB implies  IORB is singly occupied
! ICONF(I) = -IORB  implies that IORB is doubly occupied
!
! Find lexical address
!
! IF IDOREO /= 0, IREO is used to reorder lexical number
! Jeppe Olsen, November 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: ILEX_FOR_CONF
integer(kind=iwp), intent(in) :: NOCC_ORB, ICONF(NOCC_ORB), NORB, NEL, IARCW(NORB,NEL,2), IDOREO, IREO(*)
integer(kind=iwp) :: IEL, ILEX, IOCC

IEL = 0
ILEX = 1

do IOCC=1,NOCC_ORB
  if (ICONF(IOCC) > 0) then
    IEL = IEL+1
    ILEX = ILEX+IARCW(ICONF(IOCC),IEL,1)
  else if (ICONF(IOCC) < 0) then
    IEL = IEL+2
    ILEX = ILEX+IARCW(-ICONF(IOCC),IEL,2)
  end if
end do

if (IDOREO /= 0) then
  ILEX_FOR_CONF = IREO(ILEX)
else
  ILEX_FOR_CONF = ILEX
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Configuration'
call IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
write(u6,*) ' Lexical number = ',ILEX
if (IDOREO /= 0) write(u6,*) ' Reordered number = ',ILEX_FOR_CONF
#endif

end function ILEX_FOR_CONF
