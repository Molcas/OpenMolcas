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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

subroutine ClsSew()
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Real_Spherical, only: Sphere_Free
use EFP_module, only: ABC, EFP_COORS, FRAG_TYPE, lEFP
use External_Centers, only: external_centers_free
use Basis_Info, only: Basis_Info_Free, Seward_Activated
use Center_Info, only: Center_Info_Free
use Symmetry_Info, only: Symmetry_Info_Free
use SOAO_Info, only: SOAO_Info_Free
use stdalloc, only: mma_deallocate

implicit none

if (Seward_Activated) then

  call Term_Ints()
  call Free_RctFld()
  call Free_HerRW()
  call Sphere_Free()
  call SOAO_Info_Free()
  call Basis_Info_Free()
  call SYmmetry_Info_Free()
  call Center_Info_Free()
  call External_Centers_Free()
  call Free_iSD()
  call CloseR()

  if (lEFP) then
    call mma_deallocate(FRAG_TYPE)
    call mma_deallocate(ABC)
    call mma_deallocate(EFP_COORS)
    lEFP = .false.
  end if

  Seward_Activated = .false.

end if

return

end subroutine ClsSew
