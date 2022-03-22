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

subroutine TWOCT(C,S,W,THET,ENP,EPP,ABIJ,AIBJ,AJBI,BUFAB,A,B,F,FSEC,FIJKL,BUFIJ,BMN,IBMN,AC1,AC2,BUFAC)

use cpf_global, only: ICPF, IFIRST, ILIM, INCPF, INDX, IRC, ISAB, ITER, ISDCI, JSY
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: C(*), S(*), W(*), EPP(*), ABIJ(*), AIBJ(*), AJBI(*), FSEC(*)
real(kind=wp), intent(_OUT_) :: BUFAB(*), A(*), B(*), F(*), FIJKL(*), BUFIJ(*), BMN(*), AC1(*), AC2(*), BUFAC(*)
real(kind=wp), intent(in) :: THET(*), ENP(*)
integer(kind=iwp), intent(_OUT_) :: IBMN(*)

if ((ISDCI /= 0) .or. (ICPF /= 0) .or. (INCPF /= 0)) then
  ! CPF, ACPF AND SDCI
  if (ITER /= 1) then
    call DIAGC_CPF(JSY,C,S)
    if (IFIRST == 0) call ABCI_CPF(JSY,INDX,C,S,BMN,IBMN,AC1,AC2,BUFAC)
    call IJKL_CPF(JSY,INDX,C,S,FIJKL,BUFIJ,ENP,EPP)
    if (IFIRST == 0) call ABCD_CPF(JSY,INDX,ISAB,C,S,AC1,AC2,BUFAC)
  end if
  call FAIBJ_CPF(JSY,INDX,C,S,ABIJ,AIBJ,AJBI,BUFAB,A,B,F,FSEC,ENP,EPP)
else
  ! MCPF
  if (ITER /= 1) then
    call MDIAGC(JSY,C,S,W,THET,ENP,IRC(ILIM))
    if (IFIRST == 0) call MABCI(JSY,INDX,C,S,BMN,IBMN,AC1,AC2,BUFAC,W,THET,ENP,IRC(ILIM))
    call MIJKL(JSY,INDX,C,S,FIJKL,BUFIJ,W,THET,ENP,EPP,IRC(ILIM))
    if (IFIRST == 0) call MABCD(JSY,INDX,ISAB,C,S,AC1,AC2,BUFAC,W,THET,ENP,IRC(ILIM))
  end if
  call MFAIBJ(JSY,INDX,C,S,ABIJ,AIBJ,AJBI,BUFAB,A,B,F,FSEC,W,THET,ENP,EPP,IRC(ILIM))
end if

return

end subroutine TWOCT
