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

subroutine prsq(idbg,label,a,n)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idbg, n
character(len=8), intent(in) :: label
real(kind=wp), intent(in) :: a(n,n)
integer(kind=iwp) :: i, j

write(idbg,1001) label
write(idbg,1002) (j,j=1,n)
do i=1,n
  write(idbg,1003) i,(a(i,j),j=1,n)
end do

return

1001 format(' MATRIX PRINTED:',2X,A8)
1002 format(' ',4X,4(6X,I4,6X),/)
1003 format(' ',I4,4d16.8)

end subroutine prsq
