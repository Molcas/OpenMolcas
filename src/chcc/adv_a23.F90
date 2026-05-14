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

subroutine AdV_A23(V1,A,dimij,no)
! this routine does:
! A(ij,u,v) <<- V1(j,iu,v) + V1(i,jv,u)
!
! Velmi odflaknute, da sa to urobit podstatne lepsie, ale
! o4 proces osrat fok

use Index_Functions, only: iTri
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimij, no
real(kind=wp), intent(in) :: V1(no,dimij,no)
real(kind=wp), intent(inout) :: A(dimij,no,no)
integer(kind=iwp) :: i, ij, iu, j, jv, u, v

do v=1,no
  do u=1,no

    ij = 0

    do i=1,no
      iu = iTri(i,u)

      do j=1,i
        ij = ij+1
        jv = iTri(j,v)

        A(ij,u,v) = A(ij,u,v)+V1(j,iu,v)+V1(i,jv,u)

      end do
    end do

  end do
end do

return

end subroutine AdV_A23
