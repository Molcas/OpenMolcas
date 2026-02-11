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

module pt2_guga

use gugx, only: MXLEV
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: IADR10(64,2), MXCI, NG1, NG2, NG3, NG3TOT
real(kind=wp) :: CITHR, ETA(MXLEV)
character(len=8) :: CLAB10(64)

public :: CITHR, CLAB10, ETA, IADR10, MXCI, NG1, NG2, NG3, NG3TOT

end module pt2_guga
