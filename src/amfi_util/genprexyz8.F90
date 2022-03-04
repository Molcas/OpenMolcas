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

subroutine genprexyz8(preXZ)

use Constants, only: One
use Definitions, only: wp

implicit none
#include "para.fh"
real(kind=wp) :: preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)

!bs ####################################################################
!bs   additional (-) signs from the (-i) factors  in the
!bs   (-) linear combinations   (see tosigX(Y,Z).f)
!bs ####################################################################
!bs   - + - -   =>   minus
preXZ(:-1,0:,:-1,:-1) = -preXZ(:-1,0:,:-1,:-1)

return

end subroutine genprexyz8
