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

!IFG trivial
subroutine reserv_cvb(need,recn)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: need
real(kind=wp) :: recn
#include "main_cvb.fh"

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(need)
  call Unused_real(recn)
end if

end subroutine reserv_cvb
