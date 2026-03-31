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
! Copyright (C) 2022, Andy Kaiser                                      *
!***********************************************************************

module frenkel_global_vars

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: ityp, jtyp, valst
logical(kind=iwp) :: docoul, doexch, doexcitonics, excl, labb
integer(kind=iwp), allocatable :: nestla(:), nestlb(:)
real(kind=wp), allocatable :: enucb(:), vnucb(:)

public :: docoul, doexch, doexcitonics, enucb, excl, ityp, jtyp, labb, nestla, nestlb, valst, vnucb

end module frenkel_global_vars
