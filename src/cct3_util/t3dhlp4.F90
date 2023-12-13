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

subroutine t3dhlp4(w,v,dimp,dimpqr,denijk,ec,diagp,dimdiagp,addp)
! this routine realizes following procedure
! for symp=symq=symr
!
! ec = sum(pqr) [ W(pqr) . V(pqr) / Dijkpqr ]
!
! w        - W matrix (I)
! v        - V matrix (I)
! dimp     - dimension of p index (I)
! dimpqr   - dimension of q index (I)
! denijk   - sum of i,j,k diagonal (other) F parts (I)
! ec       - energy contribution
! diagp    - vector of diagonal parts of symmetry p (I)
! dimdiagp - dimension of diagonal parts in symmetry p (I)
! addp     - additional constant to p (Nocc) (I)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimpqr, dimdiagp, addp
real(kind=wp), intent(in) :: w(dimpqr), v(dimpqr), denijk, diagp(dimdiagp)
real(kind=wp), intent(out) :: ec
integer(kind=iwp) :: p, pqr, q, r
real(kind=wp) :: denijkp, denijkpq, denijkpqr

ec = Zero

pqr = 0
do p=3,dimp
  denijkp = denijk-diagp(p+addp)
  do q=2,p-1
    denijkpq = denijkp-diagp(q+addp)
    do r=1,q-1
      denijkpqr = denijkpq-diagp(r+addp)
      pqr = pqr+1
      ec = ec+w(pqr)*v(pqr)/denijkpqr
    end do
  end do
end do

return

end subroutine t3dhlp4
