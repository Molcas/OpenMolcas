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

subroutine divthelp2(t2,dima,dimb,dimi,dimj,dpa,dpb,dpi,dpj,shift_a,shift_b)
! this routine does:
! t2(a,b,i,j) = t2(a,b,i,j)/(dpi(i)+dpj(j)-dpa(a)-dpb(b))
!
! t2      - T2 matrix (I/O)
! dima    - 1 dimension of T2 (I)
! dimb    - 2 dimension of T2 (I)
! dimi    - 3 dimension of T2 (I)
! dimj    - 4 dimension of T2 (I)
! dpa     - diagonal part of Fok corresponding to irrep of a (I)
! dpb     - diagonal part of Fok corresponding to irrep of b (I)
! dpi     - diagonal part of Fok corresponding to irrep of i (I)
! dpj     - diagonal part of Fok corresponding to irrep of j (I)
! shift_a - number off occ orbitals in spin and symmetry of a (I)
! shift_b - number off occ orbitals in spin and symmetry of b (I)
!
! N.B. Since for T1 i and a are of the same spin, there is no reason
! to specify spin of dp. It must be automatically the same spin as i and a.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimi, dimj, shift_a, shift_b
real(kind=wp), intent(inout) :: t2(dima,dimb,dimi,dimj)
real(kind=wp), intent(in) :: dpa(*), dpb(*), dpi(*), dpj(*)
integer(kind=iwp) :: a, b, i, j
real(kind=wp) :: den, denij, denijb, denj

do j=1,dimj
  denj = dpj(j)
  do i=1,dimi
    denij = denj+dpi(i)
    do b=1,dimb
      denijb = denij-dpb(shift_b+b)
      do a=1,dima
        ! t2(a,b,i,j)=t2(a,b,i,j)/(denijb-dpa(shift_a+a))

        den = denijb-dpa(shift_a+a)
        if ((abs(den) >= 1.0e-7_wp) .or. (abs(t2(a,b,i,j)) > 1.0e-10_wp)) t2(a,b,i,j) = t2(a,b,i,j)/den

      end do
    end do
  end do
end do

return

end subroutine divthelp2
