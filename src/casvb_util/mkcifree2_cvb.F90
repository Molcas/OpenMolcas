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

subroutine mkcifree2_cvb(cvb,ifxstr,tconstr)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: cvb(nvb), tconstr(nvb,nvb)
integer(kind=iwp) :: ifxstr(nfxvb)

if (strucopt) then
  nfrvb = nvb
else
  nfrvb = 0
end if
nfr = nfrorb+nfrvb

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(cvb)
  call Unused_integer_array(ifxstr)
  call Unused_real_array(tconstr)
end if

end subroutine mkcifree2_cvb
