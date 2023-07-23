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

! Stuff for Cholesky vector subtraction screening (subscr):
module ChoSubScr

use Definitions, only: wp, iwp

implicit none
private

logical(kind=iwp) :: Cho_SScreen
real(kind=wp) :: SSTau, SubScrStat(2)
character(len=3) :: SSNorm
real(kind=wp), allocatable :: DSPNm(:), DSubScr(:)

public :: Cho_SScreen, DSPNm, DSubScr, SSNorm, SSTau, SubScrStat

end module ChoSubScr
