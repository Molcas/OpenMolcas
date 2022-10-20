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

subroutine t3aphlp8(a2,b,dimp,dimq,dimqr,ns,szkey)
! this routine realizes following procedure
! for symp/=symq=symr
!
! B(p,qr)=B(p,qr)+ns*(-A2(p,r,q)+A2(p,q,r))
!
! a2    - A2 matrix (I)
! b     - B  matrix (I/O)
! dimp  - dimension of p index (I)
! dimq  - dimension of q index (I)
! dimqr - dimension of qr (I)
! ns    - signum of the permutation (+-1) (I)
! szkey - set zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning

use CCT3_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimp, dimq, dimqr, ns, szkey
real(kind=wp) :: a2(dimp,dimq,dimq), b(dimp,dimqr)
integer(kind=iwp) :: nhelp, p, q, qr, qr0, r

if (szkey == 1) then
  nhelp = dimp*dimqr
  call cct3_mv0zero(nhelp,nhelp,b)
end if

if (ns == 1) then
  ! phase +1

  do q=2,dimq
    qr0 = nshf(q)
    do r=1,q-1
      qr = qr0+r
      do p=1,dimp
        b(p,qr) = b(p,qr)+a2(p,q,r)
      end do
    end do
  end do

  do q=2,dimq
    qr0 = nshf(q)
    do r=1,q-1
      qr = qr0+r
      do p=1,dimp
        b(p,qr) = b(p,qr)-a2(p,r,q)
      end do
    end do
  end do

else
  ! phase -1

  do q=2,dimq
    qr0 = nshf(q)
    do r=1,q-1
      qr = qr0+r
      do p=1,dimp
        b(p,qr) = b(p,qr)-a2(p,q,r)
      end do
    end do
  end do

  do q=2,dimq
    qr0 = nshf(q)
    do r=1,q-1
      qr = qr0+r
      do p=1,dimp
        b(p,qr) = b(p,qr)+a2(p,r,q)
      end do
    end do
  end do

end if

return

end subroutine t3aphlp8
