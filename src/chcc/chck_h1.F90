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

subroutine Chck_H1(H1,dim_,add)
! this routine test H1(i,a") = t1o(a,i)

use chcc_global, only: no, T1c
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dim_, add
real(kind=wp), intent(in) :: H1(no,dim_)
integer(kind=iwp) :: a, bad, i, ntot
real(kind=wp) :: s

bad = 0
ntot = 0

do a=1,dim_
  do i=1,no
    s = T1c(a+add,i)
    if (abs(H1(i,a)-s) > 1.0e-10_wp) then
      bad = bad+1
      !H1(i,a) = s
    end if
    ntot = ntot+1
  end do
end do

write(u6,*) ' H1 test ',bad,ntot

return

end subroutine Chck_H1
