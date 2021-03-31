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

subroutine BitMap_Localisation(PreFix)

use Localisation_globals, only: AnaNrm, ipCMO, ipMOrig, nBas, nFro, nOrb2Loc, nSym
use Index_arrays, only: iSO2Sh
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
character(len=2), intent(in) :: PreFix
#include "WrkSpc.fh"
integer(kind=iwp) :: iOff, ipCSh, ipDen, ipDSh, ipXSh, iSym, kC, kC1, kX, kX1, lCSh, lDen, lDSh, lXSh, MxBa, MxOr, n2, nBasT, &
                     nDiff, nShell
real(kind=wp) :: ThrAO
logical(kind=iwp) :: Indexation, DoF, DoG
character(len=19), parameter :: SecNam = 'BitMap_Localisation'

nBasT = nBas(1)
do iSym=2,nSym
  nBasT = nBasT+nBas(iSym)
end do

! Allocate and define some index arrays from Seward.
! --------------------------------------------------

DoF = .false.
nDiff = 0
call IniSew(DoF,nDiff)
nShell = -1
Indexation = .true.
ThrAO = Zero
DoF = .false.
DoG = .false.
call Setup_Ints(nShell,Indexation,ThrAO,DoF,DoG)
if (nShell < 1) then
  call SysAbendMsg(SecNam,'Setup_Ints failed!','nShell < 1')
end if

! Allocate max. sym. block of density matrix
! and shell based density and CMO matrices.
! ------------------------------------------

MxBa = nBas(1)
MxOr = nOrb2Loc(1)
do iSym=2,nSym
  MxBa = max(MxBa,nBas(iSym))
  MxOr = max(MxOr,nOrb2Loc(iSym))
end do
lDen = MxBa**2
lDSh = nShell**2
lCSh = nShell*MxOr
lXSh = lCSh
call GetMem('BMpLoc','Allo','Real',ipDen,lDen)
call GetMem('Dsh','Allo','Real',ipDSh,lDSh)
call GetMem('Csh','Allo','Real',ipCSh,lCSh)
call GetMem('Xsh','Allo','Real',ipXSh,lXSh)

! Compute density matrix, Den = CC^T, and set shell based matrices.
! Generate bitmap and perform sparsity analysis.
! -----------------------------------------------------------------

kC = ipMOrig
kX = ipCMO
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  call GetDens_Localisation(Work(ipDen),Work(kC1),nBas(iSym),nOrb2Loc(iSym))
  iOff = 1
  call GetSh_Localisation(Work(ipDen),nBas(iSym),nBas(iSym),Work(ipDSh),nShell,iSO2Sh(iOff),2,AnaNrm)
  call GetSh_Localisation(Work(kC1),nBas(iSym),nOrb2Loc(iSym),Work(ipCSh),nShell,iSO2Sh(iOff),1,AnaNrm)
  kX1 = kX+nBas(iSym)*nFro(iSym)
  call GetSh_Localisation(Work(kX1),nBas(iSym),nOrb2Loc(iSym),Work(ipXSh),nShell,iSO2Sh(iOff),1,AnaNrm)
  call GenBMp_Localisation(Work(ipDSh),Work(ipCSh),Work(ipXSh),nShell,iSym,'r','r','r',PreFix)
  call Anasize_Localisation(Work(ipDSh),Work(ipCSh),Work(ipXSh),nShell,nOrb2Loc(iSym),iSym)
  n2 = nBas(iSym)**2
  kC = kC+n2
  kX = kX+n2
end do
write(u6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

! De-allocations.
! ---------------

call GetMem('Xsh','Free','Real',ipXSh,lXSh)
call GetMem('Csh','Free','Real',ipCSh,lCSh)
call GetMem('Dsh','Free','Real',ipDSh,lDSh)
call GetMem('BMpLoc','Free','Real',ipDen,lDen)
DoF = .false.
DoG = .false.
call Term_Ints(DoF,DoG)

end subroutine BitMap_Localisation
