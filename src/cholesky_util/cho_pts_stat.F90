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

subroutine Cho_PTS_Stat()

use Cholesky, only: Cho_Real_Par, IntMap, LuMap, nnShl, nSym, NumCho, NumCho_G, NumChT, NumChT_G
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iTmp

if (.not. allocated(IntMap)) then
  call mma_allocate(IntMap,nnShl,Label='IntMap')
  iTmp = 0
  call IDAFile(LuMap,2,IntMap,nnShl,iTmp)
end if

if (Cho_Real_Par) then
  call iSwap(nSym,NumCho,1,NumCho_G,1)
  iTmp = NumChT
  NumChT = NumChT_G
  call Cho_Stat()
  NumChT = iTmp
  call iSwap(nSym,NumCho,1,NumCho_G,1)
else
  call Cho_Stat()
end if

if (allocated(IntMap)) call mma_deallocate(IntMap)

end subroutine Cho_PTS_Stat
