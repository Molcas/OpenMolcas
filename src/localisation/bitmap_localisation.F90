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

use Localisation_globals, only: AnaNrm, CMO, MOrig, nBas, nFro, nOrb2Loc, nSym
use iSD_Data, only: iSO2Sh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
character(len=2), intent(in) :: PreFix
integer(kind=iwp) :: iSym, kC, kC1, MxBa, MxOr, n2, nBasT, nDiff, nShell
real(kind=wp) :: ThrAO
logical(kind=iwp) :: Indexation, DoF, DoG
real(kind=wp), allocatable :: CSh(:), Den(:), DSh(:), XSh(:)
character(len=*), parameter :: SecNam = 'BitMap_Localisation'

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
call mma_allocate(Den,MxBa**2,label='BMpLoc')
call mma_allocate(DSh,nShell**2,label='Dsh')
call mma_allocate(CSh,nShell*MxOr,label='Csh')
call mma_allocate(XSh,nShell*MxOr,label='Xsh')

! Compute density matrix, Den = CC^T, and set shell based matrices.
! Generate bitmap and perform sparsity analysis.
! -----------------------------------------------------------------

kC = 1
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  call GetDens_Localisation(Den,MOrig(kC1),nBas(iSym),nOrb2Loc(iSym))
  call GetSh_Localisation(Den,nBas(iSym),nBas(iSym),DSh,nShell,iSO2Sh,2,AnaNrm)
  call GetSh_Localisation(MOrig(kC1),nBas(iSym),nOrb2Loc(iSym),CSh,nShell,iSO2Sh,1,AnaNrm)
  call GetSh_Localisation(CMO(kC1),nBas(iSym),nOrb2Loc(iSym),XSh,nShell,iSO2Sh,1,AnaNrm)
  call GenBMp_Localisation(DSh,CSh,XSh,nShell,iSym,'r','r','r',PreFix)
  call Anasize_Localisation(DSh,CSh,XSh,nShell,nOrb2Loc(iSym),iSym)
  n2 = nBas(iSym)**2
  kC = kC+n2
end do
write(u6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

! De-allocations.
! ---------------

call mma_deallocate(Den)
call mma_deallocate(CSh)
call mma_deallocate(DSh)
call mma_deallocate(XSh)
call Term_Ints()

end subroutine BitMap_Localisation
