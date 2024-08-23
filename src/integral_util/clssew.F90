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
      SubRoutine ClsSew()
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************
      use Real_Spherical, only: Sphere_Free
      use EFP_module, only: lEFP, FRAG_TYPE, ABC, EFP_COORS
      use External_Centers, only: external_centers_free
      use Basis_Info, only: Seward_Activated, Basis_Info_Free
      use Center_Info, only: Center_Info_Free
      use Symmetry_Info, only: Symmetry_Info_Free
      use SOAO_Info, only: SOAO_Info_Free
      Implicit None
!
      If (.NOT.Seward_Activated) Return
!
      Call Term_Ints()
      Call Free_RctFld()
      Call Free_HerRW()
      Call Sphere_Free()
      Call SOAO_Info_Free()
      Call Basis_Info_Free()
      Call SYmmetry_Info_Free()
      Call Center_Info_Free()
      Call External_Centers_Free()
      Call Free_iSD()
      Call CloseR()
!
      If (lEFP) Then
         Deallocate(FRAG_TYPE)
         Deallocate(ABC)
         Deallocate(EFP_COORS)
         lEFP=.FALSE.
      End If
!
      Seward_Activated=.False.
      Return
      End SubRoutine ClsSew
