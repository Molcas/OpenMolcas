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

use ipPage, only: Diskbased
use Str_Info, only: DTOC, CNSM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use MCLR_Data, only: NDTASM
use input_mclr, only: nConf, State_Sym, nCSF

implicit none
real*8 CSF(*), SD(*)
integer is
real*8, allocatable :: CTM(:)
integer iiCOPY, iprdia, iSym, i

iiCOPY = 0
iprdia = 0
nConf = max(ncsf(is),ndtasm(iS))
isym = ieor(is-1,State_Sym-1)+1
i = 2
if (isym == 1) i = 1

if (diskbased) then
  call CSDTVC_MCLR(CSF,SD,1,DTOC,CNSM(i)%ICTS,IS,iiCOPY,IPRDIA)
else
  call mma_allocate(CTM,nConf,Label='CTM')
  CTM(:) = Zero
  CTM(1:ncsf(is)) = CSF(1:ncsf(is))

  call CSDTVC_MCLR(CTM,SD,1,DTOC,CNSM(i)%ICTS,IS,iiCOPY,IPRDIA)

  call mma_deallocate(CTM)
end if

end subroutine CSF2SD
