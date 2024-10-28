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
! Copyright (C) 1990, IBM                                              *
!***********************************************************************

subroutine mHrr(la,lb,nSize,nMem)

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: la, lb
integer(kind=iwp), intent(out) :: nSize, nMem
integer(kind=iwp) :: ia, ib, nMem1, nMem2

! First find the size of the working array.

nMem = 0
nMem2 = 0
nSize = 0
do ib=0,min(la,lb)
  nMem1 = 0
  do ia=max(la,lb),la+lb-ib
    nSize = nSize+nTri_Elem1(ib)*nTri_Elem1(ia)
    nMem1 = nMem1+nTri_Elem1(ib)*nTri_Elem1(ia)
  end do
  nMem = max(nMem,nMem1+nMem2)
  nMem2 = nMem1
  if (ib == 0) nSize = 0
end do

return

end subroutine mHrr
