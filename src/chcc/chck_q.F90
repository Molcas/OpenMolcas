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

subroutine Chck_Q(Q,dimbe,addbe,dima,adda)
! check Q(be,u,i,a)

use chcc_global, only: Jc, Kc, no
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dimbe, addbe, dima, adda
real(kind=wp), intent(in) :: Q(dimbe,no,no,dima)
integer(kind=iwp) :: a, bad, be, i, u
real(kind=wp) :: s, sj, sk

bad = 0

do a=adda+1,adda+dima
  do i=1,no
    do u=1,no
      do be=addbe+1,addbe+dimbe

        sj = Jc(be,i,u,a)
        sk = Kc(i,be,u,a)

        s = 2*sj-sk

        if (abs(Q(be-addbe,u,i,a-adda)-s) > 1.0e-10_wp) then
          bad = bad+1
          !Q(be-addbe,u,i,a-adda) = s
        end if

      end do
    end do
  end do
end do

write(u6,*) ' Chck Q :',bad

return

end subroutine Chck_Q
