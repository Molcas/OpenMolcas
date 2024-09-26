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
! iterso   - second order iteration number
! MemRsv   - memory kept unallocated in LList management
! QNRTh    - threshold for QNR/C2Diis startup
! DltNTh   - convergence threshold for Norm of delta
! DltNrm   - actual Norm of delta after QNR/C2Diis extrapolation

module InfSO

use MxDM, only: MxIter
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: IterSO = 0, IterSO_Max = 30, MemRsv = 0
real(kind=wp) :: DltNrm = Zero, DltNTh = 0.2e-4_wp, Energy(MxIter) = Zero, QNRTh = 0.075_wp

public :: DltNrm, DltNth, Energy, IterSO, IterSO_Max, MemRsv, QNRTh

end module InfSO
