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

subroutine Chck_Goo(Goo)
! check Goo (i,u)

use chcc_global, only: Gooc, Hooc, no, nv, Q1, T1c
use stdalloc, only: mma_allocate
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Goo(no,no)
integer(kind=iwp) :: a, bad, i, j, u
real(kind=wp) :: s

call mma_allocate(Gooc,no,no,label='Gooc')

bad = 0

do i=1,no
  do u=1,no

    !s = Zero
    !do j=1,no
    !  do a=1,nv
    !    do b=1,nv
    !      s = s+(Two*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(a,b,u,j)
    !    end do
    !  end do
    !end do

    s = Hooc(i,u)

    do j=1,no
      do a=1,nv
        s = s+(Two*Q1(a,j,i,u)-Q1(a,i,j,u))*T1c(a,j)
      end do
    end do

    Gooc(i,u) = s

    if (abs(Goo(i,u)-s) > 1.0e-10_wp) bad = bad+1

  end do
end do

write(u6,*) ' Goo Chck :',bad

return

end subroutine Chck_Goo
