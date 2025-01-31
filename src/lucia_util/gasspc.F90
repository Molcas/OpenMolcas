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
! Copyright (C) 1998, Jeppe Olsen                                      *
!***********************************************************************

subroutine GASSPC()
! Divide orbital spaces into
!
! Inactive spaces : Orbitals that are doubly occupied in all CI spaces
! Active orbitals : Orbitals that have variable occ in atleast some spaces.
! Secondary spaces: Orbitals that are unoccupied in all spaces
!
! I_IAD : Division based upon occupation in Compound CI spaces IGSOCC
! I_IADX: Division based upon occupation in First CI space
!
! Jeppe Olsen, Summer of 98 ( not much of an summer !)

use lucia_data, only: NGAS, IGSOCC, IGSOCCX, NGSOBT
use lucia_data, only: I_IAD, I_IADX
use lucia_data, only: NELEC

implicit none
integer NEL_MAX, NEL_REF, IGAS, NTEST

! Some dummy initializtions
NEL_MAX = 0 ! jwk-cleanup

NEL_REF = NELEC(1)+NELEC(2)

! For compound space

do IGAS=1,NGAS

  if (IGAS == 1) then
    NEL_MAX = 2*NGSOBT(IGAS)
  else
    NEL_MAX = NEL_MAX+2*NGSOBT(IGAS)
  end if

  if ((IGSOCC(IGAS,1) == NEL_MAX) .and. (IGSOCC(IGAS,2) == NEL_MAX)) then
    ! Inactive  space
    I_IAD(IGAS) = 1
  else if ((IGAS > 1) .and. (IGSOCC(IGAS-1,1) == NEL_REF)) then
    ! Delete space
    I_IAD(IGAS) = 3
  else
    ! Active space
    I_IAD(IGAS) = 2
  end if

end do

! For First CI space

do IGAS=1,NGAS

  if (IGAS == 1) then
    NEL_MAX = 2*NGSOBT(IGAS)
  else
    NEL_MAX = NEL_MAX+2*NGSOBT(IGAS)
  end if

  if ((IGSOCCX(IGAS,1,1) == NEL_MAX) .and. (IGSOCCX(IGAS,2,1) == NEL_MAX)) then
    ! Inactive  space
    I_IADX(IGAS) = 1
  else if ((IGAS > 1) .and. (IGSOCCX(IGAS-1,1,1) == NEL_REF)) then
    ! Delete space
    I_IADX(IGAS) = 3
  else
    ! Active space
    I_IADX(IGAS) = 2
  end if

end do

NTEST = 0
if (NTEST >= 100) then
  write(6,*) ' Division of orbitals according to compound CI space'
  write(6,*) ' ==================================================='
  write(6,*)
  write(6,*) ' Inactive = 1, Active = 2, Delete = 3'
  write(6,*)
  call IWRTMA(I_IAD,1,NGAS,1,NGAS)
  write(6,*)
  write(6,*) ' Division of orbitals according to first CI space'
  write(6,*) ' ================================================'
  write(6,*)
  write(6,*) ' Inactive = 1, Active = 2, Delete = 3'
  write(6,*)
  call IWRTMA(I_IADX,1,NGAS,1,NGAS)
end if

end subroutine GASSPC
