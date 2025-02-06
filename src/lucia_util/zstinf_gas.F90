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

subroutine ZSTINF_GAS(IPRNT)
! Set up common block /STINF/ from information in /STINP/
!
!=========
! Input
!=========
! Information in /CGAS/ and /GASSTR/
!
!======================
! Output ( in /STINF/ )
!======================
! ISTAC (MXPSTT,2) : string type obtained by creating (ISTAC(ITYP,2))
!                    or annihilating (ISTAC(ITYP,1)) an electron
!                    from a string of type  ITYP . A zero indicates
!                    that this mapping is not included
!                    Only strings belonging to the same
!                    Orbital group are mapped
!                    mapped
! Input
! Only the first element, i.e. ISTAC  is defined

use lucia_data, only: IBGPSTR, ISTAC, MXPSTT, NGAS, NGPSTR, NGRP
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: IPRNT
integer(kind=iwp) :: IGAS, IGRP, IIGRP, MGRP, NTEST

NTEST = 0
NTEST = max(NTEST,IPRNT)
! *****************************************************************
! Mappings between strings with the same type ISTTP index, +/- 1 el
! *****************************************************************
call ISETVC(ISTAC,0,2*MXPSTT)
do IGAS=1,NGAS
  ! groups for a given gas spaces goes with increasing number of orbitals,
  ! so the first space does not have any creation mapping
  ! and the last space does not have any annihilation mapping

  MGRP = NGPSTR(IGAS)
  do IGRP=1,MGRP
    IIGRP = IGRP+IBGPSTR(IGAS)-1
    ! Annihilation map is present : IIGRP => IIGRP - 1
    if (IGRP /= 1) ISTAC(IIGRP,1) = IIGRP-1
    ! Creation map is present : IIGRP => IIGRP + 1
    if (IGRP /= MGRP) ISTAC(IIGRP,2) = IIGRP+1
  end do
end do

if (NTEST >= 10) then
  write(u6,*) ' Type - type mapping array ISTAC'
  write(u6,*) ' ==============================='
  call IWRTMA(ISTAC,NGRP,2,MXPSTT,2)
end if

end subroutine ZSTINF_GAS
