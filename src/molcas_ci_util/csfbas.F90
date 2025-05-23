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

module csfbas

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: MAXOP_LUCIA, NAEL, NBEL
integer(kind=iwp), allocatable :: CONF(:), CTS(:)

public :: CONF, CTS, MAXOP_LUCIA, NAEL, NBEL

end module csfbas
