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

subroutine ciscale2_cvb(cvec,scale,iscf,cscf)

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
dimension cvec(*)

ivec = nint(cvec(1))
iscf = 0
cscf = zero
iformat = iform_ci(ivec)
if (iformat == 0) then
  do idet=1,ndet
    work(idet+iaddr_ci(ivec)-1) = scale*work(idet+iaddr_ci(ivec)-1)
    if (abs(work(idet+iaddr_ci(ivec)-1)) > p8) then
      iscf = idet
      cscf = work(idet+iaddr_ci(ivec)-1)
    end if
  end do
else
  write(6,*) ' Unsupported format in CISCALE2 :',iformat
  call abend_cvb()
end if

return

end subroutine ciscale2_cvb
