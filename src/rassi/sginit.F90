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

subroutine SGInit(nSym,nActEl,iSpin,SGS,CIS)

use gugx, only: CIStruct, SGStruct
use MkGUGA_mod, only: MKGUGA
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nSym, nActEl, iSpin
type(SGStruct), target :: SGS
type(CIStruct) :: CIS

SGS%nSym = nSym
SGS%iSpin = iSpin
SGS%nActEl = nActEl

call MkGuga(SGS,CIS)

! Modified Arc Weights table:
call MKMAW(SGS)

! The DAW, RAW tables are no longer needed:
call mma_deallocate(SGS%RAW)
call mma_deallocate(SGS%DAW)

end subroutine SGInit
