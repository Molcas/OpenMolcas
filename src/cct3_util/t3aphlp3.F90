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

subroutine t3aphlp3(a1,a2,a3,b,dimp,dimq,dimr,dimqr,ns,szkey)
! this routine realizes following procedure
! for symp/=symq=symr
!
! B(p,qr)=B(p,qr)+ns.(A1(qr,p)-A2(p,r,q)+A3(p,q,r))
!
! a1    - A1 matrix (I)
! a2    - A2 matrix (I)
! a3    - A3 matrix (I)
! b     - B  matrix (I/O)
! dimp  - dimension of p index (I)
! dimq  - dimension of q index (I)
! dimr  - dimension of r index (I)
! dimqr - dimension of qr (I)
! ns    - signum of the permutation (+-1) (I)
! szkey - set zero key (I)
!         = 0 no vanishing
!         = 1 set B=0 at the beginning

#include "t31.fh"
integer dimp, dimq, dimr, dimqr, ns, szkey
real*8 a1(1:dimqr,1:dimp)
real*8 a2(1:dimp,1:dimr,1:dimq)
real*8 a3(1:dimp,1:dimq,1:dimr)
real*8 b(1:dimp,1:dimqr)
! help variables
integer p, q, r, qr, qr0
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
        b(p,qr) = b(p,qr)+a3(p,q,r)
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

  do qr=1,dimqr
    do p=1,dimp
      b(p,qr) = b(p,qr)+a1(qr,p)
    end do
  end do

else
  ! phase -1

  do q=2,dimq
    qr0 = nshf(q)
    do r=1,q-1
      qr = qr0+r
      do p=1,dimp
        b(p,qr) = b(p,qr)-a3(p,q,r)
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

  do qr=1,dimqr
    do p=1,dimp
      b(p,qr) = b(p,qr)-a1(qr,p)
    end do
  end do

end if

return

end subroutine t3aphlp3
