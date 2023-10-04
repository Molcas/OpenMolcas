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

subroutine cidot_cvb(cvec1,cvec2,ret)
!***********************************************************************
!*                                                                     *
!*  CIDOT  := Analogous to the blas routine DDOT                       *
!*                                                                     *
!***********************************************************************

use casvb_global, only: civbvec
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cvec1(*), cvec2(*), ret
#include "main_cvb.fh"
integer(kind=iwp) :: iformat1, iformat2, ivec1, ivec2
real(kind=wp), external :: ddot_

ivec1 = nint(cvec1(1))
ivec2 = nint(cvec2(1))
iformat1 = iform_ci(ivec1)
iformat2 = iform_ci(ivec2)
if (iformat1 /= iformat2) then
  write(u6,*) ' Format discrepancy in CIDOT :',iformat1,iformat2
  call abend_cvb()
end if
if (iformat1 == 0) then
  ret = ddot_(ndet,civbvec(:,ivec1),1,civbvec(:,ivec2),1)
else
  write(u6,*) ' Unsupported format in CIDOT :',iformat1
  call abend_cvb()
end if

return

end subroutine cidot_cvb
