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

module ChoSCF

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: ALGO, NSCREEN
real(kind=wp) :: dFKmat, dmpk = 0.045_wp
logical(kind=iwp) :: REORD

public :: ALGO, dFKmat, dmpk, NSCREEN, REORD

end module ChoSCF
