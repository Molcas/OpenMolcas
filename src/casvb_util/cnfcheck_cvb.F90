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
!****************************************************
!** Routines that deal with spatial configurations **
!****************************************************

!***********************************************************************
!*                                                                     *
!*  CNFCHECK := Check input configurations and convert to occ no       *
!*              representation if necessary.                           *
!*              Also sort in order of increasing ionicity.             *
!*                                                                     *
!***********************************************************************
subroutine cnfcheck_cvb(iconfs,nconf1,nel1)

use casvb_global, only: noe
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: nconf1, iconfs(noe,nconf1)
integer(kind=iwp), intent(in) :: nel1

call cnfcheck2_cvb(iconfs,nconf1,nel1)
call cnfsort_cvb(iconfs,nconf1,nel1)

return

end subroutine cnfcheck_cvb
