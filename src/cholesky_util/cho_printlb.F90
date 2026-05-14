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

subroutine Cho_PrintLB()

use Para_Info, only: MyRank, nProcs
use Cholesky, only: LuPri, nnBstRT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i
integer(kind=iwp), allocatable :: LB(:)

call mma_allocate(LB,[0,nProcs-1],Label='LB')
LB(:) = 0

LB(myRank) = nnBstRT(1)
call Cho_GAIGop(LB,nProcs,'+')
call Cho_Head('Cholesky vector dimension on each node','=',80,LuPri)
do i=0,nProcs-1
  write(LuPri,'(2X,A,I4,5X,A,I7)') 'Node:',i,'Dimension:',LB(i)
end do

call mma_deallocate(LB)

end subroutine Cho_PrintLB
