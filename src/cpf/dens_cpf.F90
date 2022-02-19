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

subroutine DENS_CPF(C,D,ICASE,A)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
real(kind=wp) :: C(*), D(*), A
integer(kind=iwp) :: ICASE(*)
#include "cpfmcpf.fh"
integer(kind=iwp) :: I, II, II1, ILIM, JOJ
real(kind=wp) :: EMA
integer(kind=iwp), external :: ICUNP
real(kind=r8), external :: DDOT_
! Statement function
!PAM97      EXTERNAL UNPACK
!PAM97      INTEGER UNPACK
!RL   JO(L)=IAND(ISHFT(QOCC((L+29)/30),-2*((L+29)/30*30-L)),3)
!PAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
integer(kind=iwp) :: JO, L
JO(L) = ICUNP(ICASE,L)

ILIM = NORBT*(NORBT+1)/2
call SETZ(D,ILIM)
C(IREF0) = Zero
!RL call DOTPR(C,1,C,1,A,NCONF)
A = DDOT_(NCONF,C,1,C,1)
write(u6,20) A
call XFLUSH(u6)
20 format(5X,'SUM OF SQUARED CPX(BAR)',F10.4)
C(IREF0) = One
EMA = One-A
II1 = (IREF0-1)*LN
do I=1,LN
  JOJ = JO(II1+I)
  if (JOJ >= 2) JOJ = JOJ-1
  II = I*(I+1)/2
  D(II) = JOJ*EMA
end do

return

end subroutine DENS_CPF
