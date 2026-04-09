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

subroutine DetCtl()

use MCLR_Data, only: FnCSF2SD, i12, ICISTR, idc, IDIAG, iST, LuCSF2SD, MS2, NOCSF, NOPART, pINT1, pINT2, PSSIGN
use input_mclr, only: iSpin, nActEl, nElec3, nHole1, nIrrep, nRs1, nRs2, nRs3, nSym, State_Sym
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iTmp, MNRS10, MXR4TP, MXRS30, nTRas1, nTRas3

call mma_Allocate(pINT1,nSym,Label='pInt1')
pInt1(:) = 0
call mma_Allocate(pINT2,nSym**3,Label='pInt2')
pInt2(:) = 0

NOCSF = 0
nopart = 0
nIrrep = nSym
mxr4tp = 0
idiag = 1
icistr = 1
ist = 1
i12 = 2
MS2 = iSpin-1
if (ms2 /= 0) then
  idc = 1
  pssign = Zero
else
  itmp = (ispin-1)/2
  pssign = (-One)**itmp
  idc = 2
end if

ntRas1 = sum(nRs1(1:nSym))
ntRas3 = sum(nRs3(1:nSym))
MNRS10 = max(0,2*ntRas1-nHole1)
MXRS30 = max(0,min(2*ntRas3,nElec3))
! From shells to orbitals
call ORBINF_MCLR(nSym,nRs1,nRs2,nRs3,mxr4tp) ! OK
! Number of string types
call STRTYP(ms2,nActEl,MNRS10,MXRS30)  ! looks alright
! Internal string information
call STRINF()  ! looks alright, no!
! Internal subspaces
call ICISPC(MNRS10,MXRS30)  ! looks alright
call ICISPS()  ! looks alright
! CSF information
call DANAME(LUCSF2SD,FNCSF2SD)
call CSFINF(State_sym,iSpin,1,nsym)

end subroutine DetCtl
