!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine Sort_genano(eval,evec,n,nb)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, nb
real(kind=wp), intent(inout) :: eval(n), evec(nb,n)
integer(kind=iwp) :: i, j, k, l
real(kind=wp) :: swap

do i=1,n-1
  k = i
  do j=i+1,n
   if (eval(j) > eval(k)) k = j
  end do
  if (k /= i) then
    swap    = eval(k)
    eval(k) = eval(i)
    eval(i) = swap
    do l=1,nb
      swap      =  evec(l,k)
      evec(l,k) = -evec(l,i)
      evec(l,i) =  swap
    end do
  end if
end do

return

end subroutine Sort_genano
