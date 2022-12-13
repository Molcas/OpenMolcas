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

subroutine t3aphlp4(a,b,dimp,dimpq,dimpqr,ns,szkey)
! this routine realizes following procedure
! for symp=symq=symr
!
! B(pqr)=B(pqr)+ns.(A1(qr,p)-A2(pr,q)+A3(pq,r))
!
! a      - A  matrix (I)
! b      - B  matrix (I/O)
! dimp   - dimension of p (q,r) index (I)
! dimpq  - dimension of pq index (I)
! dimpqr - dimension of pqr (I)
! ns     - signum of the permutation (+-1) (I)
! szkey  - set zero key (I)
!          = 0 no vanishing
!          = 1 set B=0 at the beginning
!
! N.B. tie sumacie by sa mozno dali trochu vylepsit

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimpq, dimpqr, ns, szkey
real(kind=wp), intent(in) :: a(dimpq,dimp)
real(kind=wp), intent(inout) :: b(dimpqr)
integer(kind=iwp) :: p, pq, pq0, pqr, pr, q, qr

if (szkey == 1) b(:) = Zero

if (ns == 1) then
  ! phase +1

  pqr = 0

  do p=3,dimp
    pq0 = (p-1)*(p-2)/2
    qr = 0
    do q=2,p-1
      pq = pq0+q
      pr = (p-1)*(p-2)/2
      b(pqr+1:pqr+q-1) = b(pqr+1:pqr+q-1)+a(qr+1:qr+q-1,p)-a(pr+1:pr+q-1,q)+a(pq,1:q-1)
      qr = qr+q-1
      pqr = pqr+q-1
    end do
  end do

else
  ! phase -1

  pqr = 0

  do p=3,dimp
    pq0 = (p-1)*(p-2)/2
    qr = 0
    do q=2,p-1
      pq = pq0+q
      pr = (p-1)*(p-2)/2
      b(pqr+1:pqr+q-1) = b(pqr+1:pqr+q-1)-a(qr+1:qr+q-1,p)+a(pr+1:pr+q-1,q)-a(pq,1:q-1)
      qr = qr+q-1
      pqr = pqr+q-1
    end do
  end do

end if

return

end subroutine t3aphlp4
