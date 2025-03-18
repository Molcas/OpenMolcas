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
function IWEYLF(NOPEN,MULTS)
! NUMBER OF CSF'S WITH NOPEN ORBITALS AND TOTAL MULTIPLICITY
! MULTS ACCORDING TO WEYLS FORMULAE
!
! (2S+1)/(NOPEN+1) * BION(NOPEN+1/0.5*NOPEN-S)

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IWEYLF
integer(kind=iwp), intent(in) :: NOPEN, MULTS
integer(kind=iwp) :: NCSF
integer(kind=iwp), external :: IBINOM

if ((NOPEN == 0) .and. (MULTS == 1)) then
  NCSF = 1
else if (mod(MULTS-1,2) /= mod(NOPEN,2)) then
  NCSF = 0
else if (mod(MULTS-1,2) == mod(NOPEN,2)) then
  NCSF = MULTS*IBINOM(NOPEN+1,(NOPEN+1-MULTS)/2)/(NOPEN+1)
end if

IWEYLF = NCSF

#ifdef _DEBUGPRINT_
write(u6,'(A,4I4)') '  IWEYLF SAYS : NOPEN MULTS NCSF : ',NOPEN,MULTS,NCSF
#endif

end function IWEYLF
