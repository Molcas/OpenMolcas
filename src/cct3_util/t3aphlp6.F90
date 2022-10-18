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

subroutine t3aphlp6(a1,a2,b,dimp,dimq,dimr,dimpq,ns,szkey)
! this routine realizes following procedure
! for symp=symq/=symr
!
! B(pq,r)=B(pq,r)+ns*(A1(q,r,p)-A2(p,r,q))
!
! a1    - A1 matrix (I)
! a2    - A2 matrix (I)
! b     - B  matrix (I/O)
! dimp  - dimension of p index (I)
! dimq  - dimension of q index (I)
! dimr  - dimension of r index (I)
! dimpq - dimension of pq (I)
! ns    - signum of the permutation (+-1) (I)
! szkey - set zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning

#include "t31.fh"
integer dimp, dimq, dimr, dimpq, ns, szkey
real*8 a1(1:dimq,1:dimr,1:dimp)
real*8 a2(1:dimp,1:dimr,1:dimq)
real*8 b(1:dimpq,1:dimr)
! help variables
integer p, q, r, pq0
integer nhelp

if (szkey == 1) then
  nhelp = dimpq*dimr
  call cct3_mv0zero(nhelp,nhelp,b)
end if

if (ns == 1) then
  ! phase +1

  do r=1,dimr
    do p=2,dimp
      pq0 = nshf(p)
      do q=1,p-1
        b(pq0+q,r) = b(pq0+q,r)-a2(p,r,q)
      end do
    end do
  end do

  do r=1,dimr
    do p=2,dimp
      pq0 = nshf(p)
      do q=1,p-1
        b(pq0+q,r) = b(pq0+q,r)+a1(q,r,p)
      end do
    end do
  end do

else
  ! phase -1

  do r=1,dimr
    do p=2,dimp
      pq0 = nshf(p)
      do q=1,p-1
        b(pq0+q,r) = b(pq0+q,r)+a2(p,r,q)
      end do
    end do
  end do

  do r=1,dimr
    do p=2,dimp
      pq0 = nshf(p)
      do q=1,p-1
        b(pq0+q,r) = b(pq0+q,r)-a1(q,r,p)
      end do
    end do
  end do

end if

return

end subroutine t3aphlp6
