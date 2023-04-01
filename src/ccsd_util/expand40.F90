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

subroutine expand40(a,b,dimpq,dimrs,dimp,dimr)
! expand a(pq,rs) -> b(p,q,r,s)
! assumption: p>q,r>s, + antisymmetry
! RISC version

integer dimpq, dimrs, dimp, dimr
real*8 a(1:dimpq,1:dimrs)
real*8 b(1:dimp,1:dimp,1:dimr,1:dimr)
! help variables
integer p, q, r, s, pq, rs
real*8 scalar

if ((dimp > 1) .and. (dimr > 1)) then

  rs = 0
  do r=2,dimr
    do s=1,r-1
      rs = rs+1

      pq = 0
      do p=2,dimp
        do q=1,p-1
          pq = pq+1

          scalar = a(pq,rs)
          b(p,q,r,s) = scalar
          b(p,q,s,r) = -scalar
          b(q,p,r,s) = -scalar
          b(q,p,s,r) = scalar

        end do
      end do
    end do
  end do

end if

do r=1,dimr
  do p=1,dimp
    do q=1,dimp
      b(p,q,r,r) = 0.0d0
    end do
  end do
end do

do s=1,dimr
  do r=1,dimr
    do p=1,dimp
      b(p,p,r,s) = 0.0d0
    end do
  end do
end do

return

end subroutine expand40
