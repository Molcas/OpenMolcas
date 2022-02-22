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

use cpf_global, only: IPRINT, IRC
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: EPP(*), EPB(*), TPQ(*), AP(*), ENP(*), T1(*), T2(*)
integer(kind=iwp) :: ICASE(*)
integer(kind=iwp) :: I, IP

IP = IRC(4)
do I=1,IP
  call TPQSET(ICASE,TPQ,I)
  call VAM(EPP,1,EPB,1,TPQ,1,T1,1,IP)
  call VDIV(ENP,1,T1,1,T2,1,IP)
  call VECSUM_CPFMCPF(T2,AP(I),IP)
  AP(I) = AP(I)*ENP(I)
end do

if (IPRINT > 5) write(u6,999) (AP(I),I=1,IP)

return

999 format(6X,'AP ',5F10.6)

end subroutine APPRIM
