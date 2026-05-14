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

!************************************
!** Print and save final VB energy **
!************************************
subroutine finalresult_cvb()

use casvb_global, only: evb, svb, variat
use Definitions, only: iwp

implicit none
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

if ((.not. variat) .and. up2date_cvb('SVB')) call add_info('SVB',[abs(svb)],1,7)
if ((.not. variat) .and. up2date_cvb('EVB')) call add_info('EVB',[evb],1,7)

return

end subroutine finalresult_cvb
