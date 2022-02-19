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

subroutine APPRIM(EPP,EPB,TPQ,AP,ENP,T1,T2,ICASE)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "cpfmcpf.fh"
dimension EPP(*), EPB(*), TPQ(*), AP(*), ENP(*), T1(*), T2(*)
dimension ICASE(*)

IP = IRC(4)
do I=1,IP
  call TPQSET(ICASE,TPQ,I)
  call VAM(EPP,1,EPB,1,TPQ,1,T1,1,IP)
  call VDIV(ENP,1,T1,1,T2,1,IP)
  call VECSUM_CPFMCPF(T2,AP(I),IP)
  AP(I) = AP(I)*ENP(I)
end do

if (IPRINT > 5) write(6,999) (AP(I),I=1,IP)
999 format(6X,'AP ',5F10.6)

return

end subroutine APPRIM
