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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimij, no
real(kind=wp) :: V1(no,dimij,no), A(dimij,no,no)
integer(kind=iwp) :: i, ij, iu, j, jv, u, v

do v=1,no
  do u=1,no

    ij = 0

    do i=1,no
      if (i >= u) then
        iu = i*(i-1)/2+u
      else
        iu = u*(u-1)/2+i
      end if

      do j=1,i
        ij = ij+1
        if (j >= v) then
          jv = j*(j-1)/2+v
        else
          jv = v*(v-1)/2+j
        end if

        A(ij,u,v) = A(ij,u,v)+V1(j,iu,v)+V1(i,jv,u)

      end do
    end do

  end do
end do

return

end subroutine AdV_A23
