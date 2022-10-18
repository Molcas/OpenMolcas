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

subroutine exth3(a,b,dimp,dimq,dimr,q,nfact)
! this routine extracts A(p,q,r) -> B_q(p,r)
!
! a     - matrix a (Input)
! b     - matrix b (Output)
! dimp  - dimension of p (Input)
! dimq  - dimension of q (Input)
! dimr  - dimension of r (Input)
! q     - value of index q (Input)
! nfact - sign (+-1,0) (Input)

integer dimp, dimq, dimr, q, nfact
real*8 a(1:dimp,1:dimq,1:dimr)
real*8 b(1:dimp,1:dimr)
! help variables
integer p, r

if (nfact == 1) then
  do r=1,dimr
    do p=1,dimp
      b(p,r) = a(p,q,r)
    end do
  end do
else if (nfact == -1) then
  do r=1,dimr
    do p=1,dimp
      b(p,r) = -a(p,q,r)
    end do
  end do
else if (nfact == 0) then
  do r=1,dimr
    do p=1,dimp
      b(p,r) = 0.0d0
    end do
  end do
end if

return

end subroutine exth3
