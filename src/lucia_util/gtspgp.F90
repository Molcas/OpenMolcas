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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

subroutine GTSPGP(IEL,ISPGP,IWAY)
! Relation between number of electrons in AS1, AS2 ... and
! supergoup number
!
! IWAY = 1 :
! Get ISPGP : Supergroup of strings
!             with IEL(*)  electrons in each AS
! IWAY = 2 :
! GET IEL(*)  : Number of electrons in each AS for supergroup ISPGP
!
!
! Jeppe Olsen, Another lonely night in Lund
!               GAS version July 1995

use lucia_data, only: NELFSPGP, NGAS, NTSPGP
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: IEL(NGAS), ISPGP
integer(kind=iwp), intent(in) :: IWAY
integer(kind=iwp) :: IGAS, JSPGP, NTEST

if (IWAY == 1) then
  ! Occupation => Number
  do JSPGP=1,NTSPGP
    if (all(NELFSPGP(1:NGAS,JSPGP) /= IEL(:))) then
      ISPGP = JSPGP
      exit
    end if
  end do
else if (IWAY == 2) then
  ! Number => Occupation
  IEL(:) = NELFSPGP(1:NGAS,ISPGP)
end if

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Output from GTSPGP'
  write(u6,*) ' IWAY ISPGP IEL ',IWAY,ISPGP,(IEL(IGAS),IGAS=1,NGAS)
end if

end subroutine GTSPGP
