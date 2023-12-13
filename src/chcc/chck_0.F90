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

subroutine Chck_0(dim_,A)
! check zero

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: dim_
real(kind=wp), intent(in) :: A(dim_)
integer(kind=iwp) :: bad, i

bad = 0
do i=1,dim_
  if (abs(A(i)) > 1.0e-10_wp) bad = bad+1
end do

write(u6,*) ' Nonzero elements ',bad,dim_

return

end subroutine Chck_0
