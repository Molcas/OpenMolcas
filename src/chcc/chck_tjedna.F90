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

subroutine Chck_Tjedna(T1)
! check T1

use chcc_global, only: no, nv, T1c
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: T1(nv,no)
integer(kind=iwp) :: a, bad, u
real(kind=wp) :: s

bad = 0

do u=1,no
  do a=1,nv
    s = T1c(a,u)
    if (abs(T1(a,u)-s) > 1.0e-10_wp) then
      bad = bad+1
      T1(a,u) = s
    end if
  end do
end do

write(u6,*) ' Tjedna   Chck :',bad

return

end subroutine Chck_Tjedna
