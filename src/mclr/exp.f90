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
Module Exp
      Implicit None
      Logical :: NewPre=.True.
      Integer ipvt,iphx,iplst
      Integer :: nexp=0, nexp_max=100
      Real*8, Allocatable:: H0S(:)
      Integer, Allocatable:: H0F(:), SBIDT(:)

      Contains
      Subroutine Exp_Close()
#include "stdalloc.fh"
      If (Allocated(H0S)) Call mma_deallocate(H0S)
      If (Allocated(H0F)) Call mma_deallocate(H0F)
      If (Allocated(SBIDT)) Call mma_deallocate(SBIDT)
      End Subroutine Exp_Close
End Module Exp
