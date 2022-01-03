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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Chk_Unitary(irc,U,n,Thr)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: check that U is unitary.

use Constants, only: Zero, One
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp) :: irc, n
real(kind=wp) :: U(n,n), Thr
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ip0, ipUTU, lUTU, n2
real(kind=wp) :: RMS, x2
real(kind=r8), external :: ddot_

if (n < 1) then
  irc = 0
  return
end if

n2 = n**2
lUTU = n2
call GetMem('UTU','Allo','Real',ipUTU,lUTU)

call dCopy_(n2,[Zero],0,Work(ipUTU),1)
ip0 = ipUTU-1
do i=1,n
  Work(ip0+n*(i-1)+i) = One
end do
call DGEMM_('T','N',n,n,n,-One,U,n,U,n,One,Work(ipUTU),n)

x2 = real(n2,kind=wp)
RMS = sqrt(dDot_(n2,Work(ipUTU),1,Work(ipUTU),1)/x2)
if (RMS > Thr) then
  irc = 1
else
  irc = 0
end if

call GetMem('UTU','Free','Real',ipUTU,lUTU)

end subroutine Chk_Unitary
