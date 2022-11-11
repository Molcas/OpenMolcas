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

module nq_MO

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nMOs = 0
real(kind=wp), allocatable :: CMO(:), D1MO(:), P2_ontop(:,:), P2MO(:)

public :: CMO, D1MO, nMOs, P2_ontop, P2MO

end module nq_MO
