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

subroutine weightfl_cvb(ix,n,nel)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n, nel
integer(kind=iwp), intent(out) :: ix(0:nel,0:n)
integer(kind=iwp) :: iel, ik

ix(:,:) = 0
ix(0,0) = 1
do iel=1,nel
  do ik=max(iel-nel+n,0),min(iel,n)
    if (ik /= 0) then
      ix(iel,ik) = ix(iel-1,ik)+ix(iel-1,ik-1)
    else
      ix(iel,ik) = ix(iel-1,ik)
    end if
  end do
end do

return

end subroutine weightfl_cvb
