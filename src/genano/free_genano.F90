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

subroutine Free_GenANO()

use stdalloc, only: mma_deallocate
use Genano_globals, only: wSet, pDSym, tDsym, Ssym, Cmo, Occ, Cmo2, Occ2, Eps, BasName

implicit none

call mma_deallocate(wSet)
call mma_deallocate(pDsym)
call mma_deallocate(tDsym)
call mma_deallocate(Ssym)
call mma_deallocate(Cmo)
call mma_deallocate(Occ)
if (allocated(Cmo2)) call mma_deallocate(Cmo2)
if (allocated(Occ2)) call mma_deallocate(Occ2)
call mma_deallocate(Eps)
call mma_deallocate(BasName)

end subroutine Free_GenANO
