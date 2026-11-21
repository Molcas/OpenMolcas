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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine CSF2SD(CSF,SD,is)
! Transforms a CSF vector to slater determinants

use CI_Info, only: CI_St
use Symmetry_Info, only: Mul
use Str_Info, only: CNSM, DTOC
use MCLR_Lucia_Data, only: NDTASM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: is
real(kind=wp), intent(in) :: CSF(CI_St(iState)%nCSF(is))
real(kind=wp), intent(_OUT_) :: SD(*)
integer(kind=iwp) :: i, iiCOPY, iSym, nC, nConf
real(kind=wp), allocatable :: CTM(:)

iiCOPY = 0
nC = CI_St(iState)%nCSF(is)
nConf = max(nC,ndtasm(iS))
isym = Mul(is,CI_St(iState)%STSym)
i = 2
if (isym == 1) i = 1

call mma_allocate(CTM,nConf,Label='CTM')
CTM(1:nC) = CSF(1:nC)
CTM(nC+1:) = Zero

call CSDTVC_MCLR_1(CTM,SD,DTOC,CNSM(i)%ICTS,IS,iiCOPY)

call mma_deallocate(CTM)

end subroutine CSF2SD
