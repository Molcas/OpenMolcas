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

subroutine t3dhlp2(w,v,dimp,dimpq,dimr,denijk,ec,diagp,diagr,dimdiagp,dimdiagr,addp,addr)
! this routine realizes following procedure
! for symp=symq/=symr
!
! ec = sum(pq,r) [ W(pq,r) . V(pq,r) / Dijkpqr ]
!
! w        - W matrix (I)
! v        - V matrix (I)
! dimp     - dimension of p (q) index (I)
! dimpq    - dimension of pq index (I)
! dimr     - dimension of r index (I)
! denijk   - sum of i,j,k diagonal (other) F parts (I)
! ec       - energy contribution
! diagp    - vector of diagonal parts of symmetry p (q) (I)
! diagr    - vector of diagonal parts of symmetry r (I)
! dimdiagp - dimension of diagonal parts in symmetry p (I)
! dimdiagr - dimension of diagonal parts in symmetry r (I)
! addp     - additional constant to p (Nocc) (I)
! addr     - additional constant to r (Nocc) (I)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimpq, dimr, dimdiagp, dimdiagr, addp, addr
real(kind=wp), intent(in) :: w(dimpq,dimr), v(dimpq,dimr), denijk, diagp(dimdiagp), diagr(dimdiagr)
real(kind=wp), intent(out) :: ec
integer(kind=iwp) :: p, pq, q, r
real(kind=wp) :: denijkpqr, denijkpr, denijkr

ec = Zero

do r=1,dimr
  denijkr = denijk-diagr(addr+r)
  pq = 0
  do p=2,dimp
    denijkpr = denijkr-diagp(p+addp)
    do q=1,p-1
      pq = pq+1
      denijkpqr = denijkpr-diagp(q+addp)
      ec = ec+w(pq,r)*v(pq,r)/denijkpqr
    end do
  end do
end do

return

end subroutine t3dhlp2
