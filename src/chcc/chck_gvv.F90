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

subroutine Chck_Gvv(Gvv)
! check Gvv(be,a)

use chcc_global, only: Gvvc, Hvvc, no, nv, Q3, T1c
use stdalloc, only: mma_allocate
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Gvv(nv,nv)
integer(kind=iwp) :: a, b, bad, be, i
real(kind=wp) :: s

call mma_allocate(Gvvc,nv,nv,label='Gvvc')

bad = 0

do be=1,nv
  do a=1,nv

    !s = Zero
    !do i=1,no
    !  do j=1,no
    !    do b=1,nv
    !      s = s+(Two*Q21(a,i,b,j)-Q21(a,j,b,i))*T2c(be,b,i,j)
    !    end do
    !  end do
    !end do
    !s =-s

    s = Hvvc(be,a)

    do i=1,no
      do b=1,nv
        s = s+(Two*Q3(a,be,b,i)-Q3(b,be,a,i))*T1c(b,i)
      end do
    end do

    Gvvc(be,a) = s

    if (abs(Gvv(be,a)-s) > 1.0e-10_wp) bad = bad+1

  end do
end do

write(u6,*) ' Gvv Chck :',bad

return

end subroutine Chck_Gvv
