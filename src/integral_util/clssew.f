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
      SubRoutine ClsSew
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************
      use setup
      use Real_Spherical, only: Sphere_Free
      use EFP_module
      use External_Centers
      use Basis_Info
      use Center_Info
      Use SOAO_Info
      use Symmetry_Info, only: Symmetry_Info_Free
      use Constants
      use stdalloc
      Implicit Real*8 (A-H,O-Z)
#include "rctfld.fh"
#include "nsd.fh"
!
      If (.NOT.Seward_Activated) Return
!
      Call Term_Ints(.False.,.True.)
      Call Free_RctFld(iXPolType)
      Call Free_HerRW()
      Call Sphere_Free()
      Call SOAO_Info_Free()
      Call Basis_Info_Free()
      Call SYmmetry_Info_Free()
      Call Center_Info_Free()
      Call External_Centers_Free()
      Call Free_iSD()
      Call Freek2()
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
      End
