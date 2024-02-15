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

subroutine sort_KQ(N,ARR,rank,proj,iopt)
! iopt = 1   => sort in ascending order
! iopt = 2   => sort in descending order

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, iopt
real(kind=wp), intent(inout) :: ARR(N)
integer(kind=iwp), intent(inout) :: rank(N), proj(N)
integer(kind=iwp) :: i, ip, ir, j
real(kind=wp) :: a

if (iopt == 1) then
  do j=2,N
    a = ARR(j)
    ir = rank(j)
    ip = proj(j)

    do i=j-1,1,-1
      if (ARR(i) <= a) exit
      ARR(i+1) = ARR(i)
      rank(i+1) = rank(i)
      proj(i+1) = proj(i)
    end do
    ARR(i+1) = a
    rank(i+1) = ir
    proj(i+1) = ip
  end do
else if (iopt == 2) then
  do j=2,N
    a = ARR(j)
    ir = rank(j)
    ip = proj(j)

    do i=j-1,1,-1
      if (ARR(i) >= a) exit
      ARR(i+1) = ARR(i)
      rank(i+1) = rank(i)
      proj(i+1) = proj(i)
    end do
    ARR(i+1) = a
    rank(i+1) = ir
    proj(i+1) = ip
  end do
else
  write(u6,'(A)') 'sort_KQ error:  iopt parameter is wrong.'
  write(u6,*) 'iopt = ',iopt
  write(u6,'(A)') 'iopt = 1, sort in ascending order'
  write(u6,'(A)') 'iopt = 2, sort in descending order'
  write(u6,'(A)') 'Return, wthout sorting'
  call xFlush(u6)
  return
end if

return

end subroutine sort_KQ
