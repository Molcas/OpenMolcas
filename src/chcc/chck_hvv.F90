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

subroutine Chck_Hvv(Hvv)
! check Hoo

use chcc_global, only: Hvvc, no, nv, Q21, T1c, T2c
use stdalloc, only: mma_allocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Hvv(nv,nv)
integer(kind=iwp) :: a, b, bad, be, i, j
real(kind=wp) :: s

call mma_allocate(Hvvc,nv,nv,label='Hvvc')

bad = 0

do be=1,nv
  do a=1,nv

    s = Zero
    do i=1,no
      do j=1,no
        do b=1,nv
          s = s+(Two*Q21(a,i,b,j)-Q21(a,j,b,i))*(T2c(be,b,i,j)+T1c(be,i)*T1c(b,j))
        end do
      end do
    end do
    s = -s

    Hvvc(be,a) = s

    if (abs(Hvv(a,be)-s) > 1.0e-10_wp) then
      bad = bad+1
      !write(u6,*) Hoo(i,u),s
    end if

  end do
end do

write(u6,*) ' Hvv Chck :',bad

return

end subroutine Chck_Hvv
