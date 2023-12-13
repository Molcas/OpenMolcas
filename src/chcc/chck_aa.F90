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

subroutine Chck_AA(A)
! check T(a,b,i,j)

use Index_Functions, only: nTri_Elem
use chcc_global, only: Ac, no
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: A(nTri_Elem(no),no,no)
integer(kind=iwp) :: bad, i, ij, j, u, v
real(kind=wp) :: s

bad = 0
do v=1,no
  do u=1,no
    ij = 0
    do i=1,no
      do j=1,i
        ij = ij+1

        s = Ac(i,j,u,v)

        if (abs(A(ij,u,v)-s) > 1.0e-10_wp) then
          bad = bad+1
          !A(ij,u,v) = s
        end if

      end do
    end do
  end do
end do

write(u6,*) ' Chck AA :',bad

return

end subroutine Chck_AA
