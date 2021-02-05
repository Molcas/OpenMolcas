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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(in) :: qInt(nInter)
real(kind=wp), intent(out) :: Hessian(nInter,nInter)
real(kind=wp), allocatable :: qInt_s(:), Hessian_s(:,:)

#ifdef _DEBUGPRINT_
call RecPrt('KHL: qInt',' ',qInt,1,nInter)
#endif

call mma_allocate(qInt_s,nInter,label='qInt_s')
call mma_allocate(Hessian_s,nInter,nInter,label='Hessian_s')

call Trans_K(qInt,qInt_s,nInter,1)
call Hessian_kriging(qInt_s,Hessian,nInter)
call BackTrans_K(Hessian,Hessian_s,nInter,nInter)
call BackTrans_Kt(Hessian_s,Hessian,nInter,nInter)

call mma_deallocate(Hessian_s)
call mma_deallocate(qInt_s)

#ifdef _DEBUGPRINT_
call RecPrt('KHL: Hessian',' ',Hessian,nInter,nInter)
#endif

end subroutine Hessian_Kriging_Layer
