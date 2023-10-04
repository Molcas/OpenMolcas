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

subroutine cinorm_cvb(cvec,cnrm)
!***********************************************************************
!*                                                                     *
!*  CINORM  := Finds the sum of squared vector elements                *
!*                                                                     *
!***********************************************************************

use casvb_global, only: civbvec
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: cvec(*), cnrm
#include "main_cvb.fh"
integer(kind=iwp) :: iformat, ivec
real(kind=wp) :: ddot_

ivec = nint(cvec(1))
iformat = iform_ci(ivec)
if (iformat == 0) then
  cnrm = ddot_(ndet,civbvec(:,ivec),1,civbvec(:,ivec),1)
else
  write(u6,*) ' Unsupported format in CINORM :',iformat
  call abend_cvb()
end if

return

end subroutine cinorm_cvb
