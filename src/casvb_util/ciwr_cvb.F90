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

!***********************************************************************
!*                                                                     *
!*  CIWR   := Write CI vector.                                         *
!*                                                                     *
!***********************************************************************
subroutine ciwr_cvb(cvec,recn)

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cvec(*), recn
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iformat, ioffs, ivec

ivec = nint(cvec(1))
iformat = iform_ci(ivec)
if (iformat == 0) then
  ioffs = 0
  call wris_cvb(iform_ci(ivec),1,recn,ioffs)
  call wris_cvb(icnt_ci(ivec),1,recn,ioffs)
  call wrrs_cvb(work(iaddr_ci(ivec)),ndet,recn,ioffs)
else
  write(u6,*) ' Unsupported format in CIWR :',iformat
  call abend_cvb()
end if

return

end subroutine ciwr_cvb
