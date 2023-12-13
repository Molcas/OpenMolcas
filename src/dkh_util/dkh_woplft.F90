!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine dkh_woplft(n,ifodd,wr,rw,p1,p2,q1,q2,t1,t2)
! Product of W(nw)P(np)=Q(np+nw)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n
logical(kind=iwp), intent(in) :: ifodd
real(kind=wp), intent(in) :: wr(n,n), rw(n,n), p1(n,n), p2(n,n)
real(kind=wp), intent(inout) :: q1(n,n), q2(n,n)
real(kind=wp), intent(out) :: t1(n,n), t2(n,n)

#include "macros.fh"
unused_var(ifodd)
unused_var(p1)
unused_var(p2)

call dmxma(n,'N','N',wr,q2,t1,One)
call dmxma(n,'N','N',rw,q1,t2,One)
q1(:,:) = t1(:,:)
q2(:,:) = t2(:,:)

return

end subroutine dkh_woplft
