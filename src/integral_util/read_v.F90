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

subroutine Read_v(lunit,work,istrt,iend,inc,ierr)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lUnit, iStrt, iEnd, Inc
real(kind=wp), intent(inout) :: work(iend)
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: i, istatus

ierr = 0
read(lunit,*,iostat=istatus) (work(istrt+i),i=0,iend-istrt,inc)
if (istatus > 0) ierr = 1

return

end subroutine Read_v
