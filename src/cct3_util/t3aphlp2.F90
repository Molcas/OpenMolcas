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

subroutine t3aphlp2(a1,a2,a3,b,dimp,dimq,dimr,dimpq,ns,szkey)
! this routine realizes following procedure
! for symp=symq/=symr
!
! B(pq,r)=B(pq,r)+ns.(A1(q,r,p)-A2(p,r,q)+A3(pq,r))
!
! a1    - A1 matrix (I)
! a2    - A2 matrix (I)
! a3    - A3 matrix (I)
! b     - B  matrix (I/O)
! dimp  - dimension of p index (I)
! dimq  - dimension of q index (I)
! dimr  - dimension of r index (I)
! dimpq - dimension of pq (I)
! ns    - signum of the permutation (+-1) (I)
! szkey - set zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning

use CCT3_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimp, dimq, dimr, dimpq, ns, szkey
real(kind=wp) :: a1(dimq,dimr,dimp), a2(dimp,dimr,dimq), a3(dimpq,dimr), b(dimpq,dimr)
integer(kind=iwp) :: nhelp, p, pq0, r

if (szkey == 1) then
  nhelp = dimpq*dimr
  call cct3_mv0zero(nhelp,nhelp,b)
end if

if (ns == 1) then
  ! phase +1

  b(:,:) = b+a3

  do r=1,dimr
    do p=2,dimp
      pq0 = nshf(p)
      b(pq0+1:pq0+p-1,r) = b(pq0+1:pq0+p-1,r)-a2(p,r,1:p-1)
    end do
  end do

  do p=2,dimp
    pq0 = nshf(p)
    b(pq0+1:pq0+p-1,:) = b(pq0+1:pq0+p-1,:)+a1(1:p-1,:,p)
  end do

else
  ! phase -1

  b(:,:) = b-a3

  do r=1,dimr
    do p=2,dimp
      pq0 = nshf(p)
      b(pq0+1:pq0+p-1,r) = b(pq0+1:pq0+p-1,r)+a2(p,r,1:p-1)
    end do
  end do

  do p=2,dimp
    pq0 = nshf(p)
    b(pq0+1:pq0+p-1,:) = b(pq0+1:pq0+p-1,:)-a1(1:p-1,:,p)
  end do

end if

return

end subroutine t3aphlp2
