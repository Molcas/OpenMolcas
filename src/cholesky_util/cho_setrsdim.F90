!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CHO_SETRSDIM(NDIMRS,MSYM,MRED,IRED,ILOC)
!
! Purpose: set reduced set dimension.

implicit real*8(a-h,o-z)
integer NDIMRS(MSYM,MRED)
#include "cholesky.fh"

if (IRED <= MAXRED) call ICOPY(NSYM,NNBSTR(1,ILOC),1,NDIMRS(1,IRED),1)

end subroutine CHO_SETRSDIM
