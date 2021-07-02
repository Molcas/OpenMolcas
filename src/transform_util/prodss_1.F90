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

subroutine ProdsS_1(AB,iAB,CMO,nMO,Y)

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "stdalloc.fh"
real*8 AB(iAB*(iAB+1)/2), CMO(iAB,nMO), Y(iAB,nMO)
real*8, allocatable :: ABSq(:)

call mma_allocate(ABSq,iAB*iAB,Label='ABSq')
call SQUARE(AB,ABSq,1,iAB,iAB)
call DGEMM_('N','N',iAB,nMO,iAB,1.0d0,ABSq,iAB,CMO,iAB,0.0d0,Y,iAB)
call mma_deallocate(ABSq)

return

end subroutine ProdsS_1
