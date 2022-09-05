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
      SubRoutine Compute_V12(V,V12,nDim)
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"
      Real*8 V(nDim,nDim), V12(nDim,nDim)

      Real*8, Allocatable :: Vec(:), VTri(:)
!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
!
      Call mma_allocate(Vec,nDim**2,Label='Vec')
      Call mma_allocate(VTri,nDim*(nDim+1)/2,Label='VTri')
!
      Call Compute_V12_(V,V12,VTri,Vec,nDim)
!
      Call mma_deallocate(VTri)
      Call mma_deallocate(Vec)
!
      Return
      End
