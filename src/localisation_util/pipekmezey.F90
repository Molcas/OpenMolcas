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
! Copyright (C) Yannick Carissan                                       *
!               Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine PipekMezey(Functional,CMO,Thrs,ThrRot,ThrGrad,BName,nBas,nOrb2Loc,nFro,nSym,nAtoms,nMxIter,Maximisation,Converged, &
                      Debug,Silent)
! Author: Y. Carissan [modified by T.B. Pedersen].
!
! Purpose: Pipek-Mezey localisation of occupied orbitals.

use OneDat, only: sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
real(kind=wp), intent(out) :: Functional
real(kind=wp), intent(inout) :: CMO(*)
real(kind=wp), intent(in) :: Thrs, ThrRot, ThrGrad
character(len=LenIn8), intent(in) :: BName(*) ! dimension should be tot. #bf
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb2Loc(nSym), nFro(nSym), nAtoms, nMxIter
logical(kind=iwp), intent(in) :: Maximisation, Debug, Silent
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: iComp, iOpt, irc, iSyLbl, kOffC, lOaux, nBasT, nFroT, nOrb2LocT
integer(kind=iwp), allocatable :: nBas_per_Atom(:), nBas_Start(:)
real(kind=wp), allocatable :: Oaux(:), Ovlp(:,:), PA(:,:,:)
character(len=8) :: Label
character(len=*), parameter :: SecNam = 'PipekMezey'

! Symmetry is NOT allowed!!
! -------------------------

if (nSym /= 1) then
  call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
end if

! Initializations.
! ----------------

Functional = -huge(Functional)

nBasT = nBas(1)
nOrb2LocT = nOrb2Loc(1)
nFroT = nFro(1)

Converged = .false.

! Read overlap matrix.
! --------------------

lOaux = nBasT*(nBasT+1)/2+4
call mma_allocate(Ovlp,nBasT,nBasT,label='Ovlp')
call mma_allocate(Oaux,lOaux,label='AuxOvlp')

irc = -1
iOpt = ibset(0,sNoOri)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(irc,iOpt,Label,iComp,Oaux,iSyLbl)
if (irc /= 0) then
  write(u6,*) SecNam,': RdOne returned ',irc
  write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
end if

if (Debug) then
  write(u6,*)
  write(u6,*) ' Triangular overlap matrix at start'
  write(u6,*) ' ----------------------------------'
  call TriPrt('Overlap',' ',Oaux,nBasT)
end if

call Tri2Rec(Oaux,Ovlp,nBasT,Debug)
call mma_deallocate(Oaux)

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

call mma_allocate(nBas_per_Atom,nAtoms,label='nB_per_Atom')
call mma_allocate(nBas_Start,nAtoms,label='nB_Start')
call BasFun_Atom(nBas_per_Atom,nBas_Start,BName,nBasT,nAtoms,Debug)

! Allocate PA array.
! ------------------
call mma_Allocate(PA,nOrb2LocT,nOrb2LocT,nAtoms,Label='PA')
PA(:,:,:) = Zero

! Localise orbitals.
! ------------------

kOffC = nBasT*nFroT+1
call PipekMezey_Iter(Functional,CMO(kOffC),Ovlp,Thrs,ThrRot,ThrGrad,PA,nBas_per_Atom,nBas_Start,BName,nBasT,nOrb2LocT,nAtoms, &
                     nMxIter,Maximisation,Converged,Debug,Silent)

! De-allocations.
! ---------------

call mma_deallocate(PA)
call mma_deallocate(nBas_per_Atom)
call mma_deallocate(nBas_Start)
call mma_deallocate(Ovlp)

end subroutine PipekMezey
