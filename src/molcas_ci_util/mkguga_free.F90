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

subroutine MKGUGA_FREE(SGS,CIS,EXS)
! PURPOSE: FREE THE GUGA TABLES

use gugx, only: CIStruct, EXStruct, SGStruct
use stdalloc, only: mma_deallocate

#include "intent.fh"

implicit none
type(SGStruct), intent(_OUT_) :: SGS
type(CIStruct), intent(_OUT_) :: CIS
type(EXStruct), intent(_OUT_) :: EXS

call mma_deallocate(SGS%ISM,safe='*')
call mma_deallocate(SGS%DRT0,safe='*')
call mma_deallocate(SGS%DOWN0,safe='*')
call mma_deallocate(SGS%DRT,safe='*')
call mma_deallocate(SGS%DOWN,safe='*')
call mma_deallocate(SGS%UP,safe='*')
call mma_deallocate(SGS%MAW,safe='*')
call mma_deallocate(SGS%LTV,safe='*')
call mma_deallocate(SGS%DAW,safe='*')
call mma_deallocate(SGS%RAW,safe='*')
call mma_deallocate(SGS%SCR,safe='*')
call mma_deallocate(SGS%Ver,safe='*')
nullify(SGS%DRTP,SGS%DOWNP)

call mma_deallocate(CIS%NOW,safe='*')
call mma_deallocate(CIS%IOW,safe='*')
call mma_deallocate(CIS%NCSF,safe='*')
call mma_deallocate(CIS%NOCSF,safe='*')
call mma_deallocate(CIS%IOCSF,safe='*')
call mma_deallocate(CIS%ICase,safe='*')

call mma_deallocate(EXS%NOCP,safe='*')
call mma_deallocate(EXS%IOCP,safe='*')
call mma_deallocate(EXS%ICoup,safe='*')
call mma_deallocate(EXS%VTab,safe='*')
call mma_deallocate(EXS%SGTMP,safe='*')
call mma_deallocate(EXS%MVL,safe='*')
call mma_deallocate(EXS%MVR,safe='*')
call mma_deallocate(EXS%USGN,safe='*')
call mma_deallocate(EXS%LSGN,safe='*')

end subroutine MKGUGA_FREE
