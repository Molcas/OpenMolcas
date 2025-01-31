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

integer function IWEYLF(NOPEN,MULTS)
! NUMBER OF CSF'S WITH NOPEN ORBITALS AND TOTAL MULTIPLICITY
! MULTS ACCORDING TO WEYLS FORMULAE
!
!     (2S+1)/(NOPEN+1) * BION(NOPEN+1/0.5NOPEN-S)

implicit real*8(A-H,O-Z)

NTEST = 0

if ((NOPEN == 0) .and. (MULTS == 1)) then
  NCSF = 1
elseif (mod(MULTS-1,2) /= mod(NOPEN,2)) then
  NCSF = 0
elseif (mod(MULTS-1,2) == mod(NOPEN,2)) then
  NCSF = MULTS*IBION_LUCIA(NOPEN+1,(NOPEN+1-MULTS)/2)/(NOPEN+1)
end if

IWEYLF = NCSF

if (NTEST /= 0) write(6,'(A,4I4)') '  IWEYLF SAYS : NOPEN MULTS NCSF : ',NOPEN,MULTS,NCSF

end function IWEYLF
