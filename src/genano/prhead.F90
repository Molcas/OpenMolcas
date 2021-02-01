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

subroutine PrHead(str,iOpt)

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: str
integer(kind=iwp), intent(in) :: iOpt
integer(kind=iwp) :: i, k, n
character(len=72) :: line
character :: fill

n = len(str)
line(1:1) = '*'
line(72:72) = '*'
fill = ' '
!-----------------------------------------------------------------------
if (iand(iOpt,1) == 0) then
  fill = '*'
  n = 0
end if
!-----------------------------------------------------------------------
k = (73-n)/2
do i=2,k
  line(i:i) = fill
end do
do i=k+n+1,71
  line(i:i) = fill
end do
if (n > 0) line(k+1:k+n) = str
!-----------------------------------------------------------------------
write(u6,*) line

return

end subroutine PrHead
