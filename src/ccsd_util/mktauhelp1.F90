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

subroutine mktauhelp1(t2,t11,t12,dima,dimb,dimi,dimj,fact)
! this routine does:
! t2(a,b,i,j) = t2(a,b,i,j) + fact [T11(a,i).T12(b,j)]
!
! t2      - T2 matrix (I/O)
! t11     - T1 amplitudes corresponding to spin ia (I)
! t12     - T1 amplitudes corresponding to spin jb (I)
! dima    - 1 dimension of T2 (I)
! dimb    - 2 dimension of T2 (I)
! dimi    - 3 dimension of T2 (I)
! dimj    - 4 dimension of T2 (I)
! fact    - numerical factor (I)
!
! N.B. symi must be syma and symj must be symb

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimi, dimj
real(kind=wp), intent(inout) :: t2(dima,dimb,dimi,dimj)
real(kind=wp), intent(in) :: t11(dima,dimi), t12(dimb,dimj), fact
integer(kind=iwp) :: b, i, j

do j=1,dimj
  do i=1,dimi
    do b=1,dimb
      t2(:,b,i,j) = t2(:,b,i,j)+fact*(t11(:,i)*t12(b,j))
    end do
  end do
end do

return

end subroutine mktauhelp1
