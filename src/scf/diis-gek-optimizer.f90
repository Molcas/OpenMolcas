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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
Subroutine DIIS_GEK_Optimizer()
!***********************************************************************
!                                                                      *
!     Object: Direct-inversion-in-the-iterative-subspace gradient-     *
!             enhanced kriging optimization.                           *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry -- BMC,                 *
!             University of Uppsala, SWEDEN                            *
!             May '22                                                  *
!***********************************************************************
use InfSO , only: iterso
use InfSCF, only: iter, mOV
use LnkLst, only: SCF_V, Init_LLs, LLx, LLGrad
Implicit None
#include "stdalloc.fh"

Integer i, j, ipq, ipg
Integer, External:: LstPtr
Real*8, Allocatable:: q(:,:), g(:,:)


Write (6,*) 'Enter DIIS-GEK Optimizer'
If (.NOT.Init_LLs) Then
   Write (6,*) 'Link list not initiated'
   Call Abend()
End If

Call mma_allocate(q,mOV,iterso,Label='q')
Call mma_allocate(g,mOV,iterso,Label='g')

j = 0
Do i = iter-iterso+1, iter
   Write (6,*) 'i,iter=',i,iter
   j = j + 1

!  Coordinates
   ipq=LstPtr(i  ,LLx)
   q(:,j)=SCF_V(ipq)%A(:)

!  Gradients
   ipg=LstPtr(i  ,LLGrad)
   g(:,j)=SCF_V(ipg)%A(:)

End Do
Call RecPrt('q',' ',q,mOV,iterso)
Call RecPrt('g',' ',g,mOV,iterso)

Call mma_deallocate(q)
Call mma_deallocate(g)

Write (6,*) 'Exit DIIS-GEK Optimizer'
End Subroutine DIIS_GEK_Optimizer
