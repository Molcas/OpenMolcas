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

subroutine Mk_iSO2Ind(iSO2Sh,iSO2Ind,nSO,nShell)

#include "stdalloc.fh"
integer iSO2Sh(nSO), iSO2Ind(nSO)
integer, allocatable :: nTemp(:)

call mma_allocate(nTemp,nShell,Label='nTemp')
call Mk_iSO2Ind_Inner(iSO2Sh,iSO2Ind,nSO,nTemp,nShell)
call mma_deallocate(nTemp)

return

end subroutine Mk_iSO2Ind
