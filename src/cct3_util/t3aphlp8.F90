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

subroutine t3aphlp8(a2,a3,b,dimp,dimq,dimqr,ns,szkey)
! this routine realizes following procedure
! for symp/=symq=symr
!
! B(p,qr)=B(p,qr)+ns*(-A2(p,r,q)+A3(p,q,r))
!
! a2    - A2 matrix (I)
! a3    - A3 matrix (I)
! b     - B  matrix (I/O)
! dimp  - dimension of p index (I)
! dimq  - dimension of q index (I)
! dimqr - dimension of qr (I)
! ns    - signum of the permutation (+-1) (I)
! szkey - set zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning

#include "t31.fh"
integer dimp, dimq, dimqr, ns, szkey
real*8 a2(1:dimp,1:dimq,1:dimq)
real*8 a3(1:dimp,1:dimq,1:dimq)
real*8 b(1:dimp,1:dimqr)
! help variables
integer p, q, r, qr0, qr
integer nhelp

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
! Avoid unused argument warnings
if (.false.) call Unused_real_array(a3)

end subroutine t3aphlp8
