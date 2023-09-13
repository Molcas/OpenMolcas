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

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension vecstr(nvb), idelstr(nzrvb)

! Zero coefficients specified by idelstr:
if (lzrvb == 0) then
  do i=1,nzrvb
    if (idelstr(i) > 0) vecstr(idelstr(i)) = 0d0
  end do
else
  ! if here, idelstr specifies which structures to *keep*:
  if (nzrvb >= 1) call fzero(vecstr,idelstr(1)-1)
  do ikeep=1,nzrvb-1
    call fzero(vecstr(idelstr(ikeep)+1),idelstr(ikeep+1)-idelstr(ikeep)-1)
  end do
end if

return

end subroutine symtrizcvb3_cvb
!*****************************************
!** All <-> constrained transformations **
!*****************************************
