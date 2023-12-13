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

subroutine MkV_A1(Ve,V1,dimo2,no)
! this routine does:
! Ve(ij,u,v) <<- V1(iu|jv)

use Index_Functions, only: iTri
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimo2, no
real(kind=wp), intent(inout) :: Ve(dimo2,no,no)
real(kind=wp), intent(in) :: V1(dimo2,dimo2)
integer(kind=iwp) :: i, ij, iu, j, jv, u, v

do v=1,no
  do u=1,no

    ij = 0
    do i=1,no
      iu = iTri(i,u)
      do j=1,i
        ij = ij+1

        jv = iTri(j,v)

        Ve(ij,u,v) = Ve(ij,u,v)+V1(iu,jv)

      end do
    end do

  end do
end do

return

end subroutine MkV_A1
