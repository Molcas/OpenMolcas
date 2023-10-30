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

subroutine cinorm2_cvb(cvec,cnrm)

use casvb_global, only: iform_ci, ndet
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: cvec(0:ndet)
real(kind=wp), intent(out) :: cnrm
integer(kind=iwp) :: iformat, ivec
real(kind=wp), external :: dnrm2_

ivec = nint(cvec(0))
iformat = iform_ci(ivec)
if (iformat == 0) then
  cnrm = dnrm2_(ndet,cvec(1:),1)
else
  write(u6,*) ' Unsupported format in CINORM2 :',iformat
  call abend_cvb()
end if

return

end subroutine cinorm2_cvb
