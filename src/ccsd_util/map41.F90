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

subroutine map41(a,b,dimp,dimq,dimr,dims,p,q,r,s,nfact)
! maping A(p1,q1,r1,s1) -> nfact* B(p2,q2,r2,s2)

real*8 a(*)
real*8 b(*)
integer dim_(4)
integer dimp, dimq, dimr, dims, p, q, r, s, nfact

dim_(p) = dimp
dim_(q) = dimq
dim_(r) = dimr
dim_(s) = dims
call map42(a,b,dimp,dimq,dimr,dims,dim_(1),dim_(2),dim_(3),dim_(4),p,q,r,s,nfact)

return

end subroutine map41
