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
! Copyright (C) 1993,1995, Jeppe Olsen                                 *
!***********************************************************************

subroutine NEWTYP(INSPGP,IACOP,ITPOP,OUTSPGP)
! an input  supergroup is given.
! apply an string of elementary operators to this supergroup and
! obtain supergroup mumber of new string
!
! Jeppe Olsen, October 1993
! GAS-version : July 95
!
! Input
! -----
!
! INSPGP  : input super group (given occupation in each AS)
! IACOP = 1 : operator is an annihilation operator
!       = 2 : operator is a creation operator
! ITPOP : orbitals space of operator
!
! Output
! ------
! OUTSPGP  : supergroup of resulting string

use lucia_data, only: NGAS, SPGPAN, SPGPCR
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: INSPGP, IACOP, ITPOP
integer(kind=iwp), intent(out) :: OUTSPGP

if (IACOP == 1) then
  OUTSPGP = SPGPAN(ITPOP+NGAS*(INSPGP-1))
else
  OUTSPGP = SPGPCR(ITPOP+NGAS*(INSPGP-1))
end if

end subroutine NEWTYP
