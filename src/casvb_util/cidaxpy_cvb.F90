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
!*  CIDAXPY := Analogous to the blas routine DAXPY                     *
!*                                                                     *
!***********************************************************************
subroutine cidaxpy_cvb(fac,cvec1,cvec2)

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cvec1(*), cvec2(*)
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iformat, ivec1, ivec2
real(kind=wp) :: fac

ivec1 = nint(cvec1(1))
ivec2 = nint(cvec2(1))
iformat = iform_ci(ivec1)
if (iformat == 0) then
  call daxpy_(ndet,fac,work(iaddr_ci(ivec1)),1,work(iaddr_ci(ivec2)),1)
else
  write(u6,*) ' Unsupported format in CIDAXPY :',iformat
  call abend_cvb()
end if
call setcnt2_cvb(ivec2,0)

return

end subroutine cidaxpy_cvb
