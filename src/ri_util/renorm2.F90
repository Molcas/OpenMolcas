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
! Copyright (C) 2008, Roland Lindh                                     *
!***********************************************************************

subroutine ReNorm2(iCnttp)

use RI_glob, only: iOffA
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iCnttp
integer(kind=iwp) :: ire_do

iOffA(:,:) = 0
do ire_do=1,2

  call ReNorm2_Inner(iCnttp)

end do

return

end subroutine ReNorm2
