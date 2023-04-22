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

subroutine mktauhelp2(t2,t1,dimab,dimij,dima,dimi,fact)
! this routine does:
! t2(ab,ij) = t2(ab,ij) + fact* (T1(a,i).t1(b,j)-t1(b,i).t1(a,j))
! for spin and symmetry of all indices equal
!
! t2    - T2 matrix (I/O)
! t1    - t1 matrix (I)
! dimab - 1 dimension of T2 (I)
! dimij - 3 dimension of T2 (I)
! dima  - number of a in this spin and symm of a (I)
! dimi  - number of a in this spin and symm of i (I)
! fact  - numerical factor (I)
!
! N.B. Since for T1 i and a are of the same spin, there is no reason
! to specify spin of dp. It must be automatically the same spin as i and a.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimab, dimij, dima, dimi
real(kind=wp), intent(inout) :: t2(dimab,dimij)
real(kind=wp), intent(in) :: t1(dima,dimi), fact
integer(kind=iwp) :: a, ab, b, i, ij, j

ij = 0
do i=2,dimi
  do j=1,i-1
    ij = ij+1

    ab = 0
    do a=2,dima
      do b=1,a-1
        ab = ab+1
        t2(ab,ij) = t2(ab,ij)+fact*(t1(a,i)*t1(b,j)-t1(b,i)*t1(a,j))

      end do
    end do
  end do
end do

return

end subroutine mktauhelp2
