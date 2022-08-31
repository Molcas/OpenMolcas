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
! Copyright (C) 1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine CtrlMO(moip,nAcO)

use Symmetry_Info, only: nIrrep
use Definitions, only: Iwp

implicit none
#include "etwas.fh"
integer(kind=iwp), intent(out) :: moip(0:nIrrep-1), nAcO
integer(kind=iwp) :: iIrrep, iTot

iTot = 0
do iIrrep=0,nIrrep-1
  moip(iIrrep) = iTot
  iTot = iTot+nAsh(iIrrep)
end do
nACO = iTot

return

end subroutine CtrlMO
