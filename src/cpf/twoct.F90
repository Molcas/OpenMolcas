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

subroutine TWOCT(H)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
dimension H(*)

call TWOCT_INTERNAL(H)

! This is to allow type punning without an explicit interface
contains
subroutine TWOCT_INTERNAL(H)
  use iso_c_binding
  real*8, target :: H(*)
  integer, pointer :: iH2(:), iH3(:), iH4(:), iH39(:), iH47(:), iH51(:)
  ILIM = 4
  if (IFIRST /= 0) ILIM = 2
  if ((ISDCI == 0) .and. (ICPF == 0) .and. (INCPF == 0)) GO TO 30
  ! CPF, ACPF AND SDCI
  if (ITER == 1) GO TO 25
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call DIAGC_CPF(iH2,H(LW(26)),H(LW(27)))
  nullify(iH2)
  if (IFIRST /= 0) GO TO 15
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(51))),iH51,[1])
  call ABCI(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(50)),iH51,H(LW(52)),H(LW(53)),H(LW(54)))
  nullify(iH2,iH3,iH51)
15 call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(47))),iH47,[1])
  call IJKL_CPF(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(46)),H(LW(47)),iH47,H(LW(31)),H(LW(32)))
  nullify(iH2,iH3,iH47)
  if (IFIRST /= 0) GO TO 25
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(4))),iH4,[1])
  call ABCD(iH2,iH3,iH4,H(LW(26)),H(LW(27)),H(LW(57)),H(LW(58)),H(LW(59)))
  nullify(iH2,iH3,iH4)
25 call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(39))),iH39,[1])
  call FAIBJ_CPF(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(36)),H(LW(37)),H(LW(38)),H(LW(39)),iH39,H(LW(40)),H(LW(41)),H(LW(42)),H(LW(43)), &
                 H(LW(31)),H(LW(32)))
  nullify(iH2,iH3,iH39)
  GO TO 50
  ! MCPF
30 if (ITER == 1) GO TO 45
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call MDIAGC(iH2,H(LW(26)),H(LW(27)),H(LW(28)),H(LW(29)),H(LW(31)),IRC(ILIM))
  nullify(iH2)
  if (IFIRST /= 0) GO TO 35
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(51))),iH51,[1])
  call MABCI(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(50)),iH51,H(LW(52)),H(LW(53)),H(LW(54)),H(LW(28)),H(LW(29)),H(LW(31)),IRC(ILIM))
  nullify(iH2,iH3,iH51)
35 call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(47))),iH47,[1])
  call MIJKL(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(46)),H(LW(47)),iH47,H(LW(28)),H(LW(29)),H(LW(31)),H(LW(32)),IRC(ILIM))
  nullify(iH2,iH3,iH47)
  if (IFIRST /= 0) GO TO 45
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(4))),iH4,[1])
  call MABCD(iH2,iH3,iH4,H(LW(26)),H(LW(27)),H(LW(57)),H(LW(58)),H(LW(59)),H(LW(28)),H(LW(29)),H(LW(31)),IRC(ILIM))
  nullify(iH2,iH3,iH4)
45 call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(39))),iH39,[1])
  call MFAIBJ(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(36)),H(LW(37)),H(LW(38)),H(LW(39)),iH39,H(LW(40)),H(LW(41)),H(LW(42)),H(LW(43)), &
              H(LW(28)),H(LW(29)),H(LW(31)),H(LW(32)),IRC(ILIM))
  nullify(iH2,iH3,iH39)
50 continue
  return
end subroutine TWOCT_INTERNAL

end subroutine TWOCT
