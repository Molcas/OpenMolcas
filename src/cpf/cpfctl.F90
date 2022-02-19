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

subroutine CPFCTL(H)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: H(*)
#include "cpfmcpf.fh"
integer(kind=iwp) :: I, IEND, IP, ISTA, ITP
real(kind=wp) :: C0, DECORR, DELE, EENP, ENER, ETOTT

call CPFCTL_INTERNAL(H)

! This is to allow type punning without an explicit interface
contains

subroutine CPFCTL_INTERNAL(H)

  real(kind=wp), target :: H(*)
  integer(kind=iwp), pointer :: iH1(:), iH2(:), iH3(:)

  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call EPSBIS(iH2,iH3,H(LW(26)),H(LW(28)),H(LW(75)))
  nullify(iH2,iH3)
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call EPSPRIM(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(32)))
  nullify(iH2,iH3)
  IP = IRC(4)
  call VECSUM_CPFMCPF(H(LW(32)),ENER,IP)
  ETOTT = ENER+POTNUC
  DELE = ETOTT-ETOT
  ETOT = ETOTT
  if (ITER == 1) then
    write(u6,'(1X,A)') ' ITER      TOTAL ENERGY          CORR ENERGY           DECREASE'
    call XFLUSH(u6)
  end if
  write(u6,'(1X,I3,3(5X,F16.8))') ITER,ETOT,ENER,DELE
  call XFLUSH(u6)
  if ((abs(DELE) < ETHRE) .and. (ITPUL /= 1)) ICONV = 1
  if ((ICONV == 0) .and. (ITER /= MAXIT)) GO TO 20
  ! If more iterations should be done, goto 20.

  if (ICONV == 1) write(u6,37)
37 format(/,5X,'CALCULATION CONVERGED')
  if (ICONV == 0) write(u6,38)
38 format(/,5X,'CALCULATION NOT COMPLETELY CONVERGED')
  if (ISDCI == 1) write(u6,30) ETOT
30 format(/,5X,'FINAL CI ENERGY',6X,F17.8)
  if (ICPF == 1) write(u6,35) ETOT
35 format(/,5X,'FINAL CPF ENERGY',5X,F17.8)
  if (INCPF == 1) write(u6,39) ETOT
39 format(/,5X,'FINAL ACPF ENERGY',4X,F17.8)
  if ((ISDCI == 0) .and. (ICPF == 0) .and. (INCPF == 0)) write(u6,36) ETOT
36 format(/,5X,'FINAL MCPF ENERGY',5X,F17.8)
  write(u6,31) ENER,POTNUC
  call XFLUSH(u6)
31 format(5X,'FINAL CORRELATION ENERGY',F14.8,'  REFERENCE ENERGY',F17.8)
  if (ISDCI == 1) call Add_Info('E_SDCI',[ETOT],1,8)
  if (ICPF == 1) call Add_Info('E_CPF',[ETOT],1,8)
  if (INCPF == 1) call Add_Info('E_ACPF',[ETOT],1,8)
  if ((ISDCI == 0) .and. (ICPF == 0) .and. (INCPF == 0)) call Add_Info('E_MCPF',[ETOT],1,8)
  call XFLUSH(u6)
  if (ISDCI == 0) GO TO 21
  EENP = H(LW(31)+IRC(4)-1)
  C0 = One/sqrt(EENP)
  DECORR = ENER*(EENP-One)
  DETOT = ETOT+DECORR
  write(u6,32) DETOT
32 format(5X,'DAVIDSON CORR. ENERGY',F17.8)
  write(u6,33) DECORR,C0
33 format(5X,'DAVIDSON CORRECTION',F19.8,'  C0 = ',F12.6)
  call XFLUSH(u6)

21 continue
  if (IPRINT > 5) then
    ISTA = LW(31)
    IEND = ISTA+IRC(4)-1
    if (IPRINT > 5) write(u6,34) (H(I),I=ISTA,IEND)
34  format(/,(5X,'ENP',5F10.6))
  end if

  return

20 continue
  ! Here if ICONV == 0 and ITER /= MAXIT (More iterations to do).
  IDIIS = 0
  if (ITPUL == MAXITP) IDIIS = 1
  call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
  call APPRIM(H(LW(32)),H(LW(75)),H(LW(30)),H(LW(76)),H(LW(31)),H(LW(79)),H(LW(80)),iH1)
  nullify(iH1)
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call CUPDATE(iH2,iH3,H(LW(26)),H(LW(27)),H(LW(76)),H(LW(33)),H(LW(80)),H(LW(31)))
  nullify(iH2,iH3)
  ITP = ITPUL+1
  call DIIS_CPF(H(LW(26)),H(LW(27)),H(LW(33)),MAXIT,H(LW(77)),ITP,H(LW(78)))

  return

end subroutine CPFCTL_INTERNAL

end subroutine CPFCTL
