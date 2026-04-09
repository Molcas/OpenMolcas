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

use Symmetry_Info, only: Mul
use Str_Info, only: CNSM, DTOC
use MCLR_Data, only: NDTASM
use input_mclr, only: nConf, nCSF, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: is
real(kind=wp), intent(in) :: CSF(nCSF(is))
real(kind=wp), intent(_OUT_) :: SD(*)
integer(kind=iwp) :: i, iiCOPY, iSym
real(kind=wp), allocatable :: CTM(:)

iiCOPY = 0
nConf = max(ncsf(is),ndtasm(iS))
isym = Mul(is,State_Sym)
i = 2
if (isym == 1) i = 1

call mma_allocate(CTM,nConf,Label='CTM')
CTM(1:ncsf(is)) = CSF(1:ncsf(is))
CTM(ncsf(is)+1:) = Zero

call CSDTVC_MCLR_1(CTM,SD,DTOC,CNSM(i)%ICTS,IS,iiCOPY)

call mma_deallocate(CTM)

end subroutine CSF2SD
