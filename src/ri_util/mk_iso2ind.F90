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
      Subroutine Mk_iSO2Ind(iSO2Sh,iSO2Ind,nSO,nShell)
#include "stdalloc.fh"
      Integer iSO2Sh(nSO), iSO2Ind(nSO)
!
      Integer, Allocatable :: nTemp(:)

      Call mma_allocate(nTemp,nShell,Label='nTemp')
      Call Mk_iSO2Ind_(iSO2Sh,iSO2Ind,nSO,nTemp,nShell)
      Call mma_deallocate(nTemp)
!
      Return
      End
