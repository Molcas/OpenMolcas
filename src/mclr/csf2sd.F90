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
use Str_Info, only: DTOC, CNSM
use MCLR_Data, only: NDTASM
use input_mclr, only: nConf, State_Sym, nCSF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero

implicit none
integer is
real*8 CSF(nCSF(is)), SD(*)
real*8, allocatable :: CTM(:)
integer iiCOPY, iSym, i

iiCOPY = 0
nConf = max(ncsf(is),ndtasm(iS))
isym = Mul(is,State_Sym)
i = 2
if (isym == 1) i = 1

call mma_allocate(CTM,nConf,Label='CTM')
CTM(1:ncsf(is)) = CSF(1:ncsf(is))
CTM(ncsf(is)+1:) = Zero

call CSDTVC_MCLR(CTM,SD,1,DTOC,CNSM(i)%ICTS,IS,iiCOPY)

call mma_deallocate(CTM)

end subroutine CSF2SD
