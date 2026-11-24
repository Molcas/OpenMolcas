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

subroutine ExplH_SA_Free()
! Release the explicit Hamiltonian of the SA CI preconditioner (ExplH_SA)

use MCLR_Data, only: ExplE, ExplV, iSel, nSel
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(iSel,safe='*')
call mma_deallocate(ExplE,safe='*')
call mma_deallocate(ExplV,safe='*')
nSel = 0

end subroutine ExplH_SA_Free
