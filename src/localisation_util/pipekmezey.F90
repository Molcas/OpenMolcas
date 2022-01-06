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
#include "WrkSpc.fh"
integer(kind=iwp) :: iComp, iOpt, ip_nBas_per_Atom, ip_nBas_Start, ipOaux, ipOvlp, irc, iSyLbl, kOffC, l_nBas_per_Atom, &
                     l_nBas_Start, lOaux, lOvlp, nBasT, nFroT, nOrb2LocT
real(kind=wp), allocatable :: PA(:,:,:)
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
lOvlp = nBasT**2
call GetMem('Ovlp','Allo','Real',ipOvlp,lOvlp)
call GetMem('AuxOvlp','Allo','Real',ipOaux,lOaux)

irc = -1
iOpt = 2
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(irc,iOpt,Label,iComp,Work(ipOaux),iSyLbl)
if (irc /= 0) then
  write(u6,*) SecNam,': RdOne returned ',irc
  write(u6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
end if

if (Debug) then
  write(u6,*)
  write(u6,*) ' Triangular overlap matrix at start'
  write(u6,*) ' ----------------------------------'
  call TriPrt('Overlap',' ',Work(ipOaux),nBasT)
end if

call Tri2Rec(Work(ipOaux),Work(ipOvlp),nBasT,Debug)
call GetMem('AuxOvlp','Free','Real',ipOaux,lOaux)

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

l_nBas_per_Atom = nAtoms
l_nBas_Start = l_nBas_per_Atom
call GetMem('nB_per_Atom','Allo','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('nB_Start','Allo','Inte',ip_nBas_Start,l_nBas_Start)
call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),BName,nBasT,nAtoms,Debug)

! Allocate PA array.
! ------------------
call mma_Allocate(PA,nOrb2LocT,nOrb2LocT,nAtoms,Label='PA')
PA(:,:,:) = Zero

! Localise orbitals.
! ------------------

kOffC = nBasT*nFroT+1
call PipekMezey_Iter(Functional,CMO(kOffC),Work(ipOvlp),Thrs,ThrRot,ThrGrad,PA,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),BName, &
                     nBasT,nOrb2LocT,nAtoms,nMxIter,Maximisation,Converged,Debug,Silent)

! De-allocations.
! ---------------

call mma_deallocate(PA)
call GetMem('nB_Start','Free','Inte',ip_nBas_Start,l_nBas_Start)
call GetMem('nB_per_Atom','Free','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('Ovlp','Free','Real',ipOvlp,lOvlp)

end subroutine PipekMezey
