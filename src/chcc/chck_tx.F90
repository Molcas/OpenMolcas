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

subroutine Chck_Tx(T)
! check T(a,b,i,j)

use chcc_global, only: no, nv, T2c
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: T(nv,no,nv,no)
integer(kind=iwp) :: a, b, bad, i, j
real(kind=wp) :: s

bad = 0
do j=1,no
  do i=1,no
    do b=1,nv
      do a=1,nv

        s = T2c(a,b,i,j)

        if (abs(T(b,i,a,j)-s) > 1.0e-10_wp) bad = bad+1

      end do
    end do
  end do
end do

write(u6,*) ' Chck T2 :',bad

return

end subroutine Chck_Tx
