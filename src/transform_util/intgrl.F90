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

module Intgrl

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: IAD2M(3,36*36), LUINTMZ, NORBZ(8), NOSHZ(8), NSYMZ

public :: IAD2M, LUINTMZ, NORBZ, NOSHZ, NSYMZ

end module Intgrl
