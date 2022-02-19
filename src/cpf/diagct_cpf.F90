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

subroutine DIAGCT_CPF(H)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
dimension H(*)

call DIAGCT_INTERNAL(H)

! This is to allow type punning without an explicit interface
contains
subroutine DIAGCT_INTERNAL(H)
  use iso_c_binding
  real*8, target :: H(*)
  integer, pointer :: iH1(:), iH2(:), iH4(:), iH11(:), iH12(:), iH13(:), iH17(:), iH18(:), iH19(:), iH20(:), iH21(:), iH22(:), &
                      iH25(:)
  ILIM = 4
  if (IFIRST /= 0) ILIM = 2
  NCONF = JSC(ILIM)
  ! Initialize sorting buffer, so that automatic detection of
  ! uninitialized variables does not give false alarms.
  call DCOPY_(LW(21)-LW(20),[0.0d0],0,H(LW(20)),1)
  ! Now H(LW(20)) and up are filled with zeroes.
  ! Similar before SORTB_CPF and SORT_CPF.
  call c_f_pointer(c_loc(H(LW(4))),iH4,[1])
  call c_f_pointer(c_loc(H(LW(20))),iH20,[1])
  call c_f_pointer(c_loc(H(LW(21))),iH21,[1])
  call c_f_pointer(c_loc(H(LW(22))),iH22,[1])
  call c_f_pointer(c_loc(H(LW(25))),iH25,[1])
  call SORTA_CPF(H(LW(20)),iH20,iH21,iH22,H(LW(10)),iH4,H(LW(25)),iH25,H(LW(23)),H(LW(24)),NINTGR)
  nullify(iH4,iH20,iH21,iH22,iH25)
  if (IFIRST == 0) then
    call DCOPY_(LW(18)-LW(17),[0.0d0],0,H(LW(17)),1)
    call c_f_pointer(c_loc(H(LW(4))),iH4,[1])
    call c_f_pointer(c_loc(H(LW(17))),iH17,[1])
    call c_f_pointer(c_loc(H(LW(18))),iH18,[1])
    call c_f_pointer(c_loc(H(LW(19))),iH19,[1])
    call SORTB_CPF(H(LW(17)),iH17,iH18,iH19,H(LW(10)),H(LW(94)),H(LW(95)),iH4,H(LW(96)))
    nullify(iH4,iH17,iH18,iH19)
  end if
  call DCOPY_(LW(12)-LW(11),[0.0d0],0,H(LW(11)),1)
  call c_f_pointer(c_loc(H(LW(11))),iH11,[1])
  call c_f_pointer(c_loc(H(LW(12))),iH12,[1])
  call c_f_pointer(c_loc(H(LW(13))),iH13,[1])
  call SORT_CPF(H(LW(11)),iH11,iH12,iH13,H(LW(14)),H(LW(15)),H(LW(16)),H(LW(10)))
  nullify(iH11,iH12,iH13)
  call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call DIAG_CPF(iH1,iH2,H(LW(11)),H(LW(14)),H(LW(15)),H(LW(16)))
  nullify(iH1,iH2)
  return
end subroutine DIAGCT_INTERNAL

end subroutine DIAGCT_CPF
