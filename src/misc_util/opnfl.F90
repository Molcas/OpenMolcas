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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

subroutine OpnFl(FName,Lu,Exists)
!***********************************************************************
!                                                                      *
! Object:                                                              *
!                                                                      *
! Called from:                                                         *
!                                                                      *
! Calling    : f_Inquire                                               *
!              isFreeUnit                                              *
!              Open                                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemical Physics,                 *
!             University of Lund, SWEDEN                               *
!             February 2000                                            *
!***********************************************************************

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: FName
integer(kind=iwp), intent(inout) :: Lu
logical(kind=iwp), intent(out) :: Exists
integer(kind=iwp), external :: isFreeUnit

! Find an unused unit number

lu = isFreeUnit(lu)
Exists = .false.

! Check that file exists

!open(unit=Lu,file=FName,status='UNKNOWN',form='FORMATTED')
call F_Inquire(FName,Exists)
call Molcas_Open(Lu,FName)

return

end subroutine OpnFl
