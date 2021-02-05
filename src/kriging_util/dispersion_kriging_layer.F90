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

subroutine Dispersion_Kriging_Layer(qInt,E_Disp,nInter)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(in) :: qInt(nInter)
real(kind=wp), intent(out) :: E_Disp
real(kind=wp), allocatable :: qInt_s(:)

call mma_allocate(qInt_s,nInter,label='qInt_s')

call Trans_K(qInt,qInt_s,nInter,1)
#ifdef _DEBUGPRINT_
call RecPrt('Dispersion_Kriging_Layer: qInt',' ',qInt,nInter,1)
call RecPrt('Dispersion_Kriging_Layer: qInt_s',' ',qInt_s,nInter,1)
#endif
call Dispersion_Kriging(qInt_s,E_Disp,nInter)

call mma_deallocate(qInt_s)

end subroutine Dispersion_Kriging_Layer
