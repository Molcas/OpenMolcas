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

subroutine Cho_X_Init_Par_GenBak()

use Para_Info, only: Is_Real_Par
use Cholesky, only: InfVec, InfVec_Bak, nSym, NumCho, NumCho_Bak
use stdalloc, only: mma_allocate

implicit none

NumCho_Bak(:) = 0
if (Is_Real_Par()) then
  call mma_allocate(InfVec_Bak,size(InfVec,1),size(InfVec,2),size(InfVec,3),Label='InfVec_Bak')
  InfVec_Bak(:,:,:) = InfVec(:,:,:)
  NumCho_Bak(1:nSym) = NumCho(1:nSym)
end if

end subroutine Cho_X_Init_Par_GenBak
