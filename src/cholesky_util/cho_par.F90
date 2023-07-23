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

subroutine Cho_ParConf(Fake)
!
! Purpose: set parallel info used in Cholesky decomposition of
!          two-electron integrals. If Fake=.True., run the
!          decomposition in serial and distribute vectors over nodes
!          at the end (no distribution is actually done, of course;
!          the relevant vectors are simply kept on the relevant
!          nodes).

use Para_Info, only: Is_Real_Par
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: Fake
#include "cho_para_info.fh"

Cho_Real_Par = Is_Real_Par() .and. (.not. Fake)

end subroutine Cho_ParConf
