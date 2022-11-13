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

subroutine t3dhlp3(w,v,dimp,dimq,dimqr,denijk,ec,diagp,diagq,dimdiagp,dimdiagq,addp,addq)
! this routine realizes following procedure
! for symp/=symq=symr
!
! ec = sum(p,qr) [ W(p,qr) . V(p,qr) / Dijkpqr ]
!
! w        - W matrix (I)
! v        - V matrix (I)
! dimp     - dimension of p index (I)
! dimq     - dimension of q (r) index (I)
! dimqr    - dimension of r index (I)
! denijk   - sum of i,j,k diagonal (other) F parts (I)
! ec       - energy contribution
! diagp    - vector of diagonal parts of symmetry p (I)
! diagq    - vector of diagonal parts of symmetry q (I)
! dimdiagp - dimension of diagonal parts in symmetry p (I)
! dimdiagq - dimension of diagonal parts in symmetry q (I)
! addp     - additional constant to p (Nocc) (I)
! addq     - additional constant to q (Nocc) (I)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimq, dimqr, dimdiagp, dimdiagq, addp, addq
real(kind=wp), intent(in) :: w(dimp,dimqr), v(dimp,dimqr), denijk, diagp(dimdiagp), diagq(dimdiagq)
real(kind=wp), intent(out) :: ec
integer(kind=iwp) :: p, q, qr, r
real(kind=wp) :: denijkpqr, denijkq, denijkqr

ec = Zero

qr = 0
do q=2,dimq
  denijkq = denijk-diagq(addq+q)
  do r=1,q-1
    denijkqr = denijkq-diagq(r+addq)
    qr = qr+1
    do p=1,dimp
      denijkpqr = denijkqr-diagp(p+addp)
      ec = ec+w(p,qr)*v(p,qr)/denijkpqr
    end do
  end do
end do

return

end subroutine t3dhlp3
