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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine symtrizcvb3_cvb(vecstr,idelstr)

use casvb_global, only: lzrvb, nvb, nzrvb
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: vecstr(nvb)
integer(kind=iwp), intent(in) :: idelstr(nzrvb)
integer(kind=iwp) :: i, ikeep

! Zero coefficients specified by idelstr:
if (lzrvb == 0) then
  do i=1,nzrvb
    if (idelstr(i) > 0) vecstr(idelstr(i)) = Zero
  end do
else
  ! if here, idelstr specifies which structures to *keep*:
  if (nzrvb >= 1) vecstr(1:idelstr(1)-1) = Zero
  do ikeep=1,nzrvb-1
    vecstr(idelstr(ikeep)+1:idelstr(ikeep+1)-1) = Zero
  end do
end if

return

end subroutine symtrizcvb3_cvb
