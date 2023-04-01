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
! dimi    - 3 domension of T2 (I)
! dimj    - 4 domension of T2 (I)
! dpa     - diagonal part of Fok corresponding to irrep of a (I)
! dpb     - diagonal part of Fok corresponding to irrep of b (I)
! dpi     - diagonal part of Fok corresponding to irrep of i (I)
! dpj     - diagonal part of Fok corresponding to irrep of j (I)
! shift_a - number off occ orbitals in spin and symmetry of a (I)
! shift_b - number off occ orbitals in spin and symmetry of b (I)
!
! N.B. Since for T1 i and a are of the same spin, there is no reason
! to specify spin of dp. It must be automatically the same spin as i and a.

integer dima, dimb, dimi, dimj, shift_a, shift_b
real*8 t2(1:dima,1:dimb,1:dimi,1:dimj)
real*8 dpa(*)
real*8 dpb(*)
real*8 dpi(*)
real*8 dpj(*)
! help variables
integer i, j, a, b
real*8 den, denj, denij, denijb

do j=1,dimj
  denj = dpj(j)
  do i=1,dimi
    denij = denj+dpi(i)
    do b=1,dimb
      denijb = denij-dpb(shift_b+b)
      do a=1,dima
        ! t2(a,b,i,j)=t2(a,b,i,j)/(denijb-dpa(shift_a+a))

        den = denijb-dpa(shift_a+a)
        if (abs(den) < 1.0d-7) then
          if (abs(t2(a,b,i,j)) > 1.0d-10) then
            t2(a,b,i,j) = t2(a,b,i,j)/den
          end if
        else
          t2(a,b,i,j) = t2(a,b,i,j)/den
        end if

      end do
    end do
  end do
end do

return

end subroutine divthelp2
