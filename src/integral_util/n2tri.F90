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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

function n2Tri(lOper)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: n2Tri
integer(kind=iwp), intent(in) :: lOper
integer(kind=iwp), external :: iPntSO

n2Tri = iPntSO(nIrrep-1,nIrrep,lOper,nBas)

return

end function n2Tri
