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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimo2, no
real(kind=wp) :: Ve(dimo2,no,no), V1(dimo2,dimo2)
integer(kind=iwp) :: i, ij, iu, j, jv, u, v

do v=1,no
  do u=1,no

    ij = 0
    do i=1,no
      do j=1,i
        ij = ij+1

        if (i > u) then
          iu = (i-1)*i/2+u
        else
          iu = (u-1)*u/2+i
        end if

        if (j > v) then
          jv = (j-1)*j/2+v
        else
          jv = (v-1)*v/2+j
        end if

        Ve(ij,u,v) = Ve(ij,u,v)+V1(iu,jv)

      end do
    end do

  end do
end do

return

end subroutine MkV_A1
