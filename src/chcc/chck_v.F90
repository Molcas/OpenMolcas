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

subroutine Chck_V(VV)
! check V

use chcc_global, only: no, nv, Q22, T1c
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: VV(nv,no,no,no)
integer(kind=iwp) :: b, bad, be, j, u, v
real(kind=wp) :: s

bad = 0

do j=1,no
  do u=1,no
    do v=1,no
      do be=1,nv

        s = Zero
        do b=1,nv
          s = s+Q22(be,b,u,j)*T1c(b,v)
        end do

        if (abs(VV(be,v,u,j)-s) > 1.0e-10_wp) then
          VV(be,v,u,j) = s
          bad = bad+1
        end if

      end do
    end do
  end do
end do

write(u6,*) ' V  Chck :',bad

return

end subroutine Chck_V
