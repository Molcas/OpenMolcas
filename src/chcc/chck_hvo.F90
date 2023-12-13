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

subroutine Chck_Hvo(Hvo)
! check Hvo

use chcc_global, only: Hvoc, no, nv, Q21, T1c
use stdalloc, only: mma_allocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Hvo(nv,no)
integer(kind=iwp) :: a, b, bad, i, j, tot
real(kind=wp) :: s

call mma_allocate(Hvoc,nv,no,label='Hvoc')

bad = 0
tot = 0

do i=1,no
  do a=1,nv

    s = Zero

    do j=1,no
      do b=1,nv
        s = s+(Two*Q21(b,j,a,i)-Q21(b,i,a,j))*T1c(b,j)
      end do
    end do

    Hvoc(a,i) = s

    if (abs(Hvo(a,i)-s) > 1.0e-10_wp) bad = bad+1
    tot = tot+1

  end do
end do

write(u6,*) ' Hvo Chck :',bad

return

end subroutine Chck_Hvo
