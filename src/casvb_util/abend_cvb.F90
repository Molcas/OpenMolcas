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

!****************
!** Error exit **
!****************
subroutine abend_cvb()

use casvb_global, only: cpu0
use Definitions, only: wp, u6

implicit none
real(kind=wp), external :: tim_cvb

write(u6,'(a)') ' Error exit CASVB.'
call date2_cvb(tim_cvb(cpu0))
call abend()

end subroutine abend_cvb
