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

subroutine DENS_CPF(C,D,ICASE,AA)

use cpf_global, only: IREF0, LN, NCONF, NORBT
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: C(*)
real(kind=wp), intent(_OUT_) :: D(*)
integer(kind=iwp), intent(in) :: ICASE(*)
real(kind=wp), intent(out) :: AA
integer(kind=iwp) :: I, II, II1, ILIM, JOJ
real(kind=wp) :: EMA
integer(kind=iwp), external :: ICUNP
real(kind=wp), external :: DDOT_

ILIM = NORBT*(NORBT+1)/2
D(1:ILIM) = Zero
C(IREF0) = Zero
AA = DDOT_(NCONF,C,1,C,1)
write(u6,20) AA
C(IREF0) = One
EMA = One-AA
II1 = (IREF0-1)*LN
do I=1,LN
  JOJ = ICUNP(ICASE,II1+I)
  if (JOJ >= 2) JOJ = JOJ-1
  II = I*(I+1)/2
  D(II) = JOJ*EMA
end do

return

20 format(5X,'SUM OF SQUARED CPX(BAR)',F10.4)

end subroutine DENS_CPF
