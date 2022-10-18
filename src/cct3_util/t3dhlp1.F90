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

subroutine t3dhlp1(w,v,dimp,dimq,dimr,denijk,ec,diagp,diagq,diagr,dimdiagp,dimdiagq,dimdiagr,addp,addq,addr)
! this routine realizes following procedure
! for symp/=symq/=symr
!
! ec = sum(p,q,r) [ W(p,q,r) . V(p,q,r) / Dijkpqr ]
!
! w        - W matrix (I)
! v        - V matrix (I)
! dimp     - dimension of p index (I)
! dimq     - dimension of q index (I)
! dimr     - dimension of r index (I)
! denijk   - sum of i,j,k diagonal (other) F parts (I)
! ec       - energy contribution
! diagp    - vector of diagonal parts of symmetry p (I)
! diagq    - vector of diagonal parts of symmetry q (I)
! diagr    - vector of diagonal parts of symmetry r (I)
! dimdiagp - dimension of diagonal parts in symmetry p (I)
! dimdiagq - dimension of diagonal parts in symmetry q (I)
! dimdiagr - dimension of diagonal parts in symmetry r (I)
! addp     - additional constant to p (Nocc) (I)
! addq     - additional constant to q (Nocc) (I)
! addr     - additional constant to r (Nocc) (I)

integer dimp, dimq, dimr
integer dimdiagp, dimdiagq, dimdiagr, addp, addq, addr
real*8 ec, denijk
real*8 w(1:dimp,1:dimq,1:dimr)
real*8 v(1:dimp,1:dimq,1:dimr)
real*8 diagp(1:dimdiagp)
real*8 diagq(1:dimdiagq)
real*8 diagr(1:dimdiagr+dimr)
! help variables
integer p, q, r
real*8 denijkr, denijkqr, denijkpqr

ec = 0.0d0

do r=1,dimr
  denijkr = denijk-diagr(addr+r)
  do q=1,dimq
    denijkqr = denijkr-diagq(q+addq)
    do p=1,dimp
      denijkpqr = denijkqr-diagp(p+addp)
      ec = ec+w(p,q,r)*v(p,q,r)/denijkpqr
    end do
  end do
end do

return

end subroutine t3dhlp1
