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

subroutine Chck_Vx(V)
! check V

use chcc_global, only: Jc, Kc, no, nv
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: V(nv,no,nv,no)
integer(kind=iwp) :: a, bad, be, i, u
real(kind=wp) :: s

bad = 0
do a=1,nv
  do i=1,no
    do u=1,no
      do be=1,nv

        s = Two*Jc(be,i,u,a)-Kc(i,be,u,a)
        s = Kc(i,be,u,a)

        if (abs(V(be,u,a,i)-s) > 1.0e-10_wp) then
          bad = bad+1
          V(be,u,a,i) = s
        end if

      end do
    end do
  end do
end do

write(u6,*) ' Chck Vx :',bad

return

end subroutine Chck_Vx
