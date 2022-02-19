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

subroutine DENSCT_CPF(H,LIC0)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: H(LIC0)
integer(kind=iwp) :: LIC0
#include "cpfmcpf.fh"
real(kind=wp) :: A

call DENSCT_INTERNAL(H)

! This is to allow type punning without an explicit interface
contains

subroutine DENSCT_INTERNAL(H)

  real(kind=wp), target :: H(*)
  integer(kind=iwp), pointer :: iH1(:), iH2(:), iH3(:), iH34(:)

  call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
  call DENS_CPF(H(LW(26)),H(LW(62)),iH1,A)
  nullify(iH1)

  ! MULTIPLY C BY MP

  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call c_f_pointer(c_loc(H(LW(34))),iH34,[1])
  call NPSET(iH2,iH3,H(LW(26)),H(LW(30)),H(LW(31)),H(LW(72)),H(LW(27)),H(LW(28)),H(LW(32)),iH34)
  nullify(iH2,iH3,iH34)

  call ONECT(H)
  if (A > One) then
    write(u6,*) 'DENSCT_CPF Error: A>1.0 (See code.)'
  end if
  call NATCT(H,LIC0)

  return

end subroutine DENSCT_INTERNAL

end subroutine DENSCT_CPF
