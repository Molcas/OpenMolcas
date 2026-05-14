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

module NAC

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: NACstates(2)
real(kind=wp) :: Ediff
logical(kind=iwp) :: DoCSF, isCSF, isNAC

public :: DoCSF, Ediff, isCSF, isNAC, NACstates

end module NAC
