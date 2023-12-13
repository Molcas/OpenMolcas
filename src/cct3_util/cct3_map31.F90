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

subroutine cct3_map31(a,b,dimp,dimq,dimr,p,q,r,nfact)
! mapping A(p1,q1,r1) -> nfact*B(p2,q2,r2)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimr, p, q, r, nfact
real(kind=wp), intent(in) :: a(*)
real(kind=wp), intent(_OUT_) :: b(*)
integer(kind=iwp) :: dim_(3)

dim_(p) = dimp
dim_(q) = dimq
dim_(r) = dimr
call cct3_map32(a,b,dimp,dimq,dimr,dim_(1),dim_(2),dim_(3),p,q,nfact)

return

end subroutine cct3_map31
