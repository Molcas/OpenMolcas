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

subroutine divthelp3(t2,dimab,dimij,dpa,dpi,dima,dimi,shift)
! this routine does:
! t2(ab,ij) = t2(ab,ij)/(dp(i)+dp(j)-dp(a)-dp(b))
! for spin and symmetry of a = b
! and spin and symmetry of i = j
!
! t2    - T2 matrix (I/O)
! dimab - 1 dimension of T2 (I)
! dimij - 3 dimension of T2 (I)
! dpa   - diagonal part of Fok corresponding to irrep of a (I)
! dpi   - diagonal part of Fok corresponding to irrep of i (I)
! dima  - number of a in this spin and symm of a (I)
! dimi  - number of a in this spin and symm of i (I)
! shift - number of occ orbitals in spin and symmetry of a (I)
!
! N.B. Since for T1 i and a are of the same spin, there is no reason
! to specify spin of dp. It must be automatically the same spin as i and a.

integer dimab, dimij, dima, dimi, shift
real*8 t2(1:dimab,1:dimij)
real*8 dpa(*)
real*8 dpi(*)
! help variables
integer i, j, a, b, ij, ab
real*8 den, deni, denij, denija

ij = 0
do i=2,dimi
  deni = dpi(i)
  do j=1,i-1
    denij = deni+dpi(j)
    ij = ij+1

    ab = 0
    do a=2,dima
      denija = denij-dpa(shift+a)
      do b=1,a-1
        ab = ab+1
        ! t2(ab,ij)=t2(ab,ij)/(denija-dpa(shift+b))

        den = denija-dpa(shift+b)
        if (abs(den) < 1.0d-7) then
          if (abs(t2(ab,ij)) > 1.0d-10) then
            t2(ab,ij) = t2(ab,ij)/den
          end if
        else
          t2(ab,ij) = t2(ab,ij)/den
        end if

      end do
    end do
  end do
end do

return

end subroutine divthelp3
