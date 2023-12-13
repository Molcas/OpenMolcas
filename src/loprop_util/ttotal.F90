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

subroutine Ttotal(T1,T2,T3,T4,Ttot,Ttot_Inv,nDim)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(in) :: T1(nDim,nDim), T2(nDim,nDim), T3(nDim,nDim), T4(nDim,nDim)
real(kind=wp), intent(out) :: Ttot(nDim,nDim), Ttot_Inv(nDim,nDim)
real(kind=wp) :: DET
real(kind=wp), allocatable :: Temp(:,:), Temp2(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Ttot=T1*T2*T3*T4

!lg write(u6,*) 'Ttotal ', nDim
call mma_allocate(Temp,nDim,nDim,label='Temp')
call mma_allocate(Temp2,nDim,nDim,label='Temp2')
call DGEMM_('N','N',nDim,nDim,nDim,One,T1,nDim,T2,nDim,Zero,Temp,nDim)
call DGEMM_('N','N',nDim,nDim,nDim,One,Temp,nDim,T3,nDim,Zero,Temp2,nDim)
call DGEMM_('N','N',nDim,nDim,nDim,One,Temp2,nDim,T4,nDim,Zero,Ttot,nDim)
call mma_deallocate(Temp)
call mma_deallocate(Temp2)
!lg call RecPrt('T_TOT',' ',Ttot,nDim,nDim)
call MINV(Ttot,Ttot_Inv,DET,nDim)

return

end subroutine Ttotal
