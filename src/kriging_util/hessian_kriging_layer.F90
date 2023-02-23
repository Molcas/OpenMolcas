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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

subroutine Hessian_Kriging_Layer(qInt,Hessian,nInter)

use kriging_mod, only: nSet
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(in) :: qInt(nInter)
real(kind=wp), intent(out) :: Hessian(nInter,nInter,nSet)
integer(kind=iwp) :: iSet
real(kind=wp), allocatable :: qInt_s(:), Hessian_s(:,:)

#ifdef _DEBUGPRINT_
call RecPrt('KHL: qInt',' ',qInt,1,nInter)
#endif

call mma_allocate(qInt_s,nInter,label='qInt_s')
call mma_allocate(Hessian_s,nInter,nInter,label='Hessian_s')

call Trans_K(qInt,qInt_s,nInter,1)
call Hessian_kriging(qInt_s,Hessian,nInter)

do iSet=1,nSet
  call BackTrans_K(Hessian(:,:,iSet),Hessian_s,nInter,nInter)
  call BackTrans_Kt(Hessian_s,Hessian(:,:,iSet),nInter,nInter)
end do

call mma_deallocate(Hessian_s)
call mma_deallocate(qInt_s)

#ifdef _DEBUGPRINT_
do iSet=1,nSet
  call RecPrt('KHL: Hessian',' ',Hessian(:,:,iSet),nInter,nInter)
end do
#endif

end subroutine Hessian_Kriging_Layer
