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

subroutine Chck_K(K,dimbe,addbe,dima,adda)
! check K(be,u,i,a)

use chcc_global, only: Kc, no
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dimbe, addbe, dima, adda
real(kind=wp), intent(in) :: K(dimbe,no,no,dima)
integer(kind=iwp) :: a, bad, be, i, u
real(kind=wp) :: s

bad = 0
do a=adda+1,adda+dima
  do i=1,no
    do u=1,no
      do be=addbe+1,addbe+dimbe

        s = Kc(i,be,u,a)

        if (abs(K(be-addbe,u,i,a-adda)-s) > 1.0e-10_wp) bad = bad+1

      end do
    end do
  end do
end do

write(u6,*) ' Chck K :',bad
!99 format(a9,1x,i8,1x,4(i3,1x))

return

end subroutine Chck_K
