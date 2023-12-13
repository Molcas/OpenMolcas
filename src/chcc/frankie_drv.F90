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

subroutine frankie_drv(NChHere)

use chcc_global, only: nfr, no, nv, printkey
use Cholesky, only: NumCho
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: NChHere

! ----------------------------------------------------------------------

! - transform AO Cholesky vectors to MO basis locally on each node
!   with fragmented cholesky index. _AI1 -> _CDtmp1

call frankie(nfr,no,nv,printkey)

! take local # of Cholesky Vectors on this node
NChHere = NumCho(1)

return

end subroutine frankie_drv
