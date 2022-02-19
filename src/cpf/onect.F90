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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine ONECT(H)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp

implicit none
real(kind=wp) H(*)
#include "cpfmcpf.fh"
integer(kind=iwp) :: ILIM

call ONECT_INTERNAL(H)

! This is to allow type punning without an explicit interface
contains

subroutine ONECT_INTERNAL(H)

  real(kind=wp), target :: H(*)
  integer(kind=iwp), pointer :: iH1(:), iH2(:), iH3(:), iH63(:)

  ILIM = 4
  if (IFIRST /= 0) ILIM = 2
  if ((ICPF == 0) .and. (ISDCI == 0) .and. (INCPF == 0)) GO TO 15
  ! CPF AND SDCI
  if (IDENS == 1) GO TO 10
  ! (AI/JK) INTEGRALS
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(63))),iH63,[1])
  call AI_CPF(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(62)),H(LW(63)),iH63,H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),H(LW(31)),H(LW(32)),1)
  nullify(iH2,iH3,iH63)
10 call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call FIJ(iH1,iH2,iH3,H(LW(26)),H(LW(27)),H(LW(62)),H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),H(LW(31)),H(LW(32)))
  nullify(iH1,iH2,iH3)
  GO TO 20
  ! MCPF
15 if (IDENS == 1) GO TO 5
  ! (AI/JK) INTEGRALS
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(63))),iH63,[1])
  call MAI(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(62)),H(LW(63)),iH63,H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),H(LW(28)),H(LW(29)), &
           H(LW(31)),H(LW(32)),IRC(ILIM),1)
  nullify(iH2,iH3,iH63)
5 call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call MFIJ(iH1,iH2,iH3,H(LW(26)),H(LW(27)),H(LW(62)),H(LW(64)),H(LW(65)),H(LW(66)),H(LW(67)),H(LW(28)),H(LW(29)),H(LW(31)), &
            H(LW(32)),IRC(ILIM))
  nullify(iH1,iH2,iH3)
20 continue

  return

end subroutine ONECT_INTERNAL

end subroutine ONECT
