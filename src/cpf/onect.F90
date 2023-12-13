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

subroutine ONECT(C,S,W,THET,ENP,EPP,FC,BUFIN,A,B,FK,DBK)

use cpf_global, only: ICASE, ICPF, IDENS, ILIM, INCPF, INDX, IRC, ISDCI, JSY
use Definitions, only: wp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: C(*), S(*), W(*), EPP(*), FC(*), FK(*)
real(kind=wp), intent(in) :: THET(*), ENP(*)
real(kind=wp), intent(_OUT_) :: BUFIN(*), A(*), B(*), DBK(*)

if ((ICPF /= 0) .or. (ISDCI /= 0) .or. (INCPF /= 0)) then
  ! CPF AND SDCI
  if (IDENS /= 1) then
    ! (AI/JK) INTEGRALS
    call AI_CPF(JSY,INDX,C,S,FC,BUFIN,A,B,FK,DBK,ENP,EPP,1)
  end if
  call FIJ_CPF(ICASE,JSY,INDX,C,S,FC,A,B,FK,DBK,ENP,EPP)
else
  ! MCPF
  if (IDENS /= 1) then
    ! (AI/JK) INTEGRALS
    call MAI(JSY,INDX,C,S,FC,BUFIN,A,B,FK,DBK,W,THET,ENP,EPP,IRC(ILIM),1)
  end if
  call MFIJ(ICASE,JSY,INDX,C,S,FC,A,B,FK,DBK,W,THET,ENP,EPP,IRC(ILIM))
end if

return

end subroutine ONECT
