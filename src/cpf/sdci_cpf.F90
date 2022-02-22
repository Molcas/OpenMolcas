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

subroutine SDCI_CPF(H,iH,LIC0)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use cpf_global, only: ICONV, ICPF, IDENS, INCPF, IPRINT, IRC, IREF0, IREST, ISDCI, ITER, ITPUL, JSC, LIC, LW, MAXIT, NREF
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: H(*)
integer(kind=iwp) :: iH(*), LIC0

call SDCI_CPF_INTERNAL(H)

! This is to allow type punning without an explicit interface
contains

subroutine SDCI_CPF_INTERNAL(H)

  real(kind=wp), target :: H(*)
  integer(kind=iwp), pointer :: iH1(:), iH2(:), iH3(:)

  LIC = LIC0
  IPRINT = 5
  IDENS = 0
  ! INPUT, SORTING AND DIAGONAL ELEMENTS
  call READIN_CPF(H,iH)
  call DIAGCT_CPF(H)
  ITER = 1
  if (IREST == 1) ITER = ITER+1
  ITPUL = 1
  if (IREST == 0) call START_CPF(H(LW(26)),JSC(4),IREF0)
  if (IREST == 1) call RESTART_CPFMCPF(H(LW(26)),JSC(4))
  if ((ICPF == 0) .and. (INCPF == 0) .and. (ISDCI == 0)) then
    call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
    call THETSET(iH1,H(LW(29)),IRC(4))
    nullify(iH1)
  end if
  do
    call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
    call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
    call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
    call NPSET(iH2,iH3,H(LW(26)),H(LW(30)),H(LW(31)),H(LW(72)),H(LW(27)),H(LW(28)),H(LW(32)),iH1)
    nullify(iH1,iH2,iH3)
    call TWOCT(H)
    call ONECT(H)
    call CPFCTL(H)
    ITER = ITER+1
    ITPUL = ITPUL+1
    if ((ITER > MAXIT) .or. (ICONV == 1)) exit
  end do
  ! FIRST ORDER DENSITY MATRIX
  IDENS = 1

  call DENSCT_CPF(H,LIC0)

  if (NREF > 1) then
    write(u6,*) ' This is a single reference program, but more than'
    write(u6,*) ' one reference state has been specified in the'
    write(u6,*) ' GUGA program. Change input to GUGA and run again.'
  end if

  return

end subroutine SDCI_CPF_INTERNAL

end subroutine SDCI_CPF
