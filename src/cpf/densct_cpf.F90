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

subroutine DENSCT_CPF(C,S,W,THET,TPQ,ENP,EPP,ICASE_,FC,BUFIN,A,B,FK,DBK,TEMP)

use cpf_global, only: ICASE, INDX, JSY
use Constants, only: One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: C(*), ENP(*), FK(*)
real(kind=wp), intent(_OUT_) :: S(*), W(*), TPQ(*), EPP(*), FC(*), BUFIN(*), A(*), B(*), DBK(*), TEMP(*)
real(kind=wp), intent(in) :: THET(*)
integer(kind=iwp), intent(in) :: ICASE_(*)
real(kind=wp) :: AA

call DENS_CPF(C,FC,ICASE,AA)

! MULTIPLY C BY MP

call NPSET(JSY,INDX,C,TPQ,ENP,TEMP,S,W,EPP,ICASE_)

call ONECT(C,S,W,THET,ENP,EPP,FC,BUFIN,A,B,FK,DBK)
if (AA > One) then
  write(u6,*) 'DENSCT_CPF Error: AA>1.0 (See code.)'
end if
call NATCT(C,FC)

return

end subroutine DENSCT_CPF
