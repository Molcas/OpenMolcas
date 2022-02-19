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

implicit real*8(A-H,O-Z)
dimension C(*), D(*)
dimension ICASE(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
! Statement function
!PAM97      EXTERNAL UNPACK
!PAM97      INTEGER UNPACK
!RL   JO(L)=IAND(ISHFT(QOCC((L+29)/30),-2*((L+29)/30*30-L)),3)
!PAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
JO(L) = ICUNP(ICASE,L)

ILIM = NORBT*(NORBT+1)/2
call SETZ(D,ILIM)
C(IREF0) = D0
!RL call DOTPR(C,1,C,1,A,NCONF)
A = DDOT_(NCONF,C,1,C,1)
write(6,20) A
call XFLUSH(6)
20 format(5X,'SUM OF SQUARED CPX(BAR)',F10.4)
C(IREF0) = D1
EMA = D1-A
II1 = (IREF0-1)*LN
do I=1,LN
  JOJ = JO(II1+I)
  if (JOJ >= 2) JOJ = JOJ-1
  II = I*(I+1)/2
  D(II) = JOJ*EMA
end do

return

end subroutine DENS_CPF
