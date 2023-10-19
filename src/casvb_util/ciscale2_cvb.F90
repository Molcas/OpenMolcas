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

subroutine ciscale2_cvb(cvec,scl,iscf,cscf)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
real(kind=wp) :: cvec(0:ndet), scl, cscf
integer(kind=iwp) :: iscf
integer(kind=iwp) :: idet, iformat, ivec

ivec = nint(cvec(0))
iscf = 0
cscf = Zero
iformat = iform_ci(ivec)
if (iformat == 0) then
  do idet=1,ndet
    cvec(idet) = scl*cvec(idet)
    if (abs(cvec(idet)) > p8) then
      iscf = idet
      cscf = cvec(idet)
    end if
  end do
else
  write(u6,*) ' Unsupported format in CISCALE2 :',iformat
  call abend_cvb()
end if

return

end subroutine ciscale2_cvb
