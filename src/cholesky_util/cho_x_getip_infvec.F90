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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Cho_X_GetIP_InfVec(InfVcT)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: get pointer to InfVec array for all vectors.

use Cholesky, only: InfVec
#ifdef _MOLCAS_MPP_
use Cholesky, only: Cho_Real_Par, InfVec_Bak
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), pointer, intent(out) :: InfVct(:,:,:)

#ifdef _MOLCAS_MPP_
if (Cho_Real_Par) then
  if (allocated(InfVec_Bak)) then
    InfVcT => InfVec_Bak
  else
    call Cho_Quit('Initialization problem in Cho_X_GetIP_InfVec',103)
  end if
else
#endif
  InfVcT => InfVec
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine Cho_X_GetIP_InfVec
