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

subroutine NEXT(P,DPS,CN)

implicit real*8(A-H,O-Z)
dimension P(*), DPS(*), CN(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"

IAD = IADDP(1)
call dDAFILE(Lu_CI,2,P,NCONF,IAD)
ITM = ITPUL-1
do I=1,ITM
  IN = I+1
  CTOT = 0.0d00
  do J=IN,ITPUL
    CTOT = CTOT+CN(J)
  end do
  IAD = IADDP(I+1)
  call dDAFILE(Lu_CI,2,DPS,NCONF,IAD)
  call VSMA(DPS,1,CTOT,P,1,P,1,NCONF)
end do
if (IPRINT >= 15) write(6,19) (P(I),I=1,NCONF)
19 format(6X,'C(NEXT)',5F10.6)

IADC(ITPUL+2) = IAD

return

end subroutine NEXT
