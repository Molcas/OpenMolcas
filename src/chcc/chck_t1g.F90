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

subroutine Chck_T1g(T1,dima,adda)
! this routine tests T1g

use chcc_global, only: no, T1c
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dima, adda
real(kind=wp), intent(inout) :: T1(no,dima)
integer(kind=iwp) :: a, bad, i, ntot
real(kind=wp) :: s

bad = 0
ntot = 0

do i=1,no
  do a=adda+1,adda+dima

    s = T1c(a,i)

    if (abs(T1(i,a-adda)-s) > 1.0e-10_wp) then
      bad = bad+1
      T1(i,a-adda) = s
    end if
    ntot = ntot+1

  end do
end do

write(u6,*) ' T1g   ',bad,ntot

return

end subroutine Chck_T1g
