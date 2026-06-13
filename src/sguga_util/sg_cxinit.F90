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

subroutine SG_CXInit(SGS,CIS,EXS)

use sguga, only: SGStruct, CIStruct, EXStruct, MkCoup, MkSeg, MkNrCoup
use stdalloc, only: mma_deallocate

implicit none
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(out) :: CIS
type(EXStruct), intent(out) :: EXS

CIS%nMidV = SGS%MVEnd-SGS%MVSta+1
CIS%nIpWlk = 1+(SGS%MidLev-1)/15
CIS%nIpWlk = max(CIS%nIpWlk,1+(SGS%nLev-SGS%MidLev-1)/15)

! Calculate segment values

call MkSeg(SGS,CIS,EXS)

call MkNrCoup(SGS,CIS,EXS)

call MkCoup(SGS,CIS,EXS)

call mma_deallocate(CIS%ISgm)
call mma_deallocate(CIS%VSgm)
call mma_deallocate(CIS%IVR)

end subroutine SG_CXInit
