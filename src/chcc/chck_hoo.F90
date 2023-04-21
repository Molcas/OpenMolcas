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

subroutine Chck_Hoo(Hoo)
! check Hoo

use chcc_global, only: Hooc, no, nv, Q21, T1c, T2c
use stdalloc, only: mma_allocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Hoo(no,no)
integer(kind=iwp) :: a, b, bad, i, j, u
real(kind=wp) :: s

call mma_allocate(Hooc,no,no,label='Hooc')

bad = 0

do i=1,no
  do u=1,no

    s = Zero
    do j=1,no
      do a=1,nv
        do b=1,nv
          s = s+(Two*Q21(a,i,b,j)-Q21(a,j,b,i))*(T2c(a,b,u,j)+T1c(a,u)*T1c(b,j))
        end do
      end do
    end do

    Hooc(i,u) = s

    if (abs(Hoo(i,u)-s) > 1.0e-10_wp) then
      bad = bad+1
      !write(u6,*) Hoo(i,u),s
    end if

  end do
end do

write(u6,*) ' Hoo Chck :',bad

return

end subroutine Chck_Hoo
