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

subroutine ncset_cvb(ic)

use casvb_global, only: istackrep
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: ic
integer(kind=iwp) :: nc_zeroed, nconvinone
logical(kind=iwp), external :: istkprobe_cvb

if (istkprobe_cvb(istackrep)) then
  call istkpop_cvb(istackrep,nc_zeroed)
  call istkpop_cvb(istackrep,nconvinone)
  if ((ic == 0) .or. (ic == 1)) then
    nconvinone = nconvinone+1
  else if (ic > 1) then
    nconvinone = 0
    nc_zeroed = 1
  else
    nconvinone = -1
    nc_zeroed = 1
  end if
  call istkpush_cvb(istackrep,nconvinone)
  call istkpush_cvb(istackrep,nc_zeroed)
end if

return

end subroutine ncset_cvb
