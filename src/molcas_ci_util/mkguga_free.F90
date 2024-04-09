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

if (allocated(SGS%ISM)) call mma_deallocate(SGS%ISM)
if (allocated(SGS%DRT0)) call mma_deallocate(SGS%DRT0)
if (allocated(SGS%DOWN0)) call mma_deallocate(SGS%DOWN0)
if (allocated(SGS%DRT)) call mma_deallocate(SGS%DRT)
if (allocated(SGS%DOWN)) call mma_deallocate(SGS%DOWN)
if (allocated(SGS%UP)) call mma_deallocate(SGS%UP)
if (allocated(SGS%MAW)) call mma_deallocate(SGS%MAW)
if (allocated(SGS%LTV)) call mma_deallocate(SGS%LTV)
if (allocated(SGS%DAW)) call mma_deallocate(SGS%DAW)
if (allocated(SGS%RAW)) call mma_deallocate(SGS%RAW)
if (allocated(SGS%SCR)) call mma_deallocate(SGS%SCR)
if (allocated(SGS%Ver)) call mma_deallocate(SGS%Ver)
SGS%DRTP => null()
SGS%DOWNP => null()

if (allocated(CIS%NOW)) call mma_deallocate(CIS%NOW)
if (allocated(CIS%IOW)) call mma_deallocate(CIS%IOW)
if (allocated(CIS%NCSF)) call mma_deallocate(CIS%NCSF)
if (allocated(CIS%NOCSF)) call mma_deallocate(CIS%NOCSF)
if (allocated(CIS%IOCSF)) call mma_deallocate(CIS%IOCSF)
if (allocated(CIS%ICase)) call mma_deallocate(CIS%ICase)

if (allocated(EXS%NOCP)) call mma_deallocate(EXS%NOCP)
if (allocated(EXS%IOCP)) call mma_deallocate(EXS%IOCP)
if (allocated(EXS%ICoup)) call mma_deallocate(EXS%ICoup)
if (allocated(EXS%VTab)) call mma_deallocate(EXS%VTab)
if (allocated(EXS%SGTMP)) call mma_deallocate(EXS%SGTMP)
if (allocated(EXS%MVL)) call mma_deallocate(EXS%MVL)
if (allocated(EXS%MVR)) call mma_deallocate(EXS%MVR)
if (allocated(EXS%USGN)) call mma_deallocate(EXS%USGN)
if (allocated(EXS%LSGN)) call mma_deallocate(EXS%LSGN)

end subroutine MKGUGA_FREE
