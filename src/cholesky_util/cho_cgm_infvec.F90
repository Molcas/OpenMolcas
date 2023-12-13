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
! Copyright (C) Thomas Bondo Pedersen                                  *
!               2020,2021, Roland Lindh                                *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Cho_CGM_InfVec(InfVcT,NVT,n)

use Definitions, only: iwp

implicit none
integer(kind=iwp), pointer, intent(out) :: InfVcT(:,:,:)
integer(kind=iwp), intent(in) :: n
integer(kind=iwp), intent(out) :: NVT(n)

call Cho_X_GetIP_InfVec(InfVcT)
call Cho_X_GetTotV(NVT,n)

end subroutine Cho_CGM_InfVec
