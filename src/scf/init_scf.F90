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

subroutine Init_SCF()

use SCF_Arrays, only: Dens, TwoHam, Vxc
use InfSCF, only: MapDns, Two_Thresholds
use RICD_Info, only: Do_DCCD
use Constants, only: Zero

implicit none
integer i, nActEl

! Clear Dens and TwoHam matrices
Dens(:,:,:) = Zero
TwoHam(:,:,:) = Zero
Vxc(:,:,:) = Zero

! Set number of active shells on the RUNFILE to zero

call Peek_iScalar('nSym',i)
nActEl = 0
call Put_iScalar('nActel',nActEl)

call IniLLs()
! clear MapDns ...
MapDns(:) = 0

Two_Thresholds = .not. Do_DCCD

return

end subroutine Init_SCF
