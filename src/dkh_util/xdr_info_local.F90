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

subroutine XDR_Info_Local(n,indx,nb,ib,map)
! calculate information of Local blocks

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n, indx(n)
integer(kind=iwp), intent(out) :: nb, ib(n), map(n)
integer(kind=iwp) :: i, j, k, m

ib(:) = 0
nb = 0
k = 0
do i=1,n
  if (ib(i) == 0) then
    nb = nb+1
    m = k
    k = k+1
    map(k) = i
    do j=i+1,n
      if ((ib(j) == 0) .and. (indx(j) == indx(i))) then
        k = k+1
        map(k) = j
        ib(j) = -1
      end if
    end do
    ib(nb) = k-m
  end if
end do
!DP write(u6,'(a,9i4)') ' nb,k,n : ',nb,k,n
!DP write(u6,'(a,9i5)') ' blocks : ',(ib(i),i=1,nb)
!DP write(u6,'(a,999i4)') ' map ; ',(map(i),i=1,n)
return

end subroutine XDR_Info_Local
