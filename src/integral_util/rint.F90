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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

function RinT(iOpT,nOpT,iOpR)
!***********************************************************************
!                                                                      *
! Object: to return .true. if R is in {T}.                             *
!                                                                      *
! Called from: TwoEl                                                   *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             June '90                                                 *
!***********************************************************************

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: RinT
integer(kind=iwp), intent(in) :: nOpT, iOpT(nOpT), iOpR
integer(kind=iwp) :: i

RinT = .false.
do i=1,nOpT
  if (iOpT(i) == iOpR) then
    RinT = .true.
    return
  end if
end do

end function RinT
