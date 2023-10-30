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

subroutine symtrizcvb_cvb(vecstr)

use casvb_global, only: iconstruc, idelstr, ipermzeta, izeta, nconstr, nvb, tconstr
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: vecstr(nvb)
real(kind=wp) :: dum(1)

if (iconstruc == 0) then
  return
else if (iconstruc == 1) then
  call symtrizcvb2_cvb(vecstr,izeta,ipermzeta)
  call symtrizcvb3_cvb(vecstr,idelstr)
else if (iconstruc == 2) then
  call schmidtd_cvb(tconstr,nconstr,vecstr,1,dum,nvb,0)
end if

return

end subroutine symtrizcvb_cvb
