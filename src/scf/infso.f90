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
! Module for second order update info
!
!     iterso   - second order iteration number
!     MemRsv   - memory kept unallocated in LList management
!     QNRTh    - threshold for QNR/C2Diis startup
!     DltNTh   - convergence threshold for Norm of delta
!     DltNrm   - actual Norm of delta after QNR/C2Diis extrapolation

Module InfSO
use MxDM, only: MxIter
use Constants, only: Zero
Private
Public :: IterSO, IterSO_Max, MemRsv, QNRTh, DltNth, DltNrm, Energy

Integer :: i

Integer :: IterSO=0
Integer :: IterSO_Max=30
Integer :: MemRsv=0
Real*8  :: QNRTh = 0.075d+00
Real*8  :: DltNTh= 0.2d-4
Real*8  :: DltNrm= Zero
Real*8  :: Energy(MxIter)=[(Zero,i=1,MxIter)]

End Module InfSO
