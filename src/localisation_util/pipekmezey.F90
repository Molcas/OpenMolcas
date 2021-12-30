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

subroutine PipekMezey(Functional,CMO,Thrs,ThrRot,ThrGrad,Name,nBas,nOrb2Loc,nFro,nSym,nAtoms,nMxIter,Maximisation,Converged,Debug, &
                      Silent)
! Author: Y. Carissan [modified by T.B. Pedersen].
!
! Purpose: Pipek-Mezey localisation of occupied orbitals.

implicit real*8(a-h,o-z)
#include "Molcas.fh"
real*8 CMO(*)
integer nBas(nSym), nOrb2Loc(nSym), nFro(nSym)
logical Maximisation, Converged, Debug, Silent
character*(LENIN8) Name(*) ! dimension should be tot. #bf
#include "WrkSpc.fh"
#include "stdalloc.fh"
real*8, allocatable :: PA(:,:,:)
character*10 SecNam
parameter(SecNam='PipekMezey')
character*8 Label

! Symmetry is NOT allowed!!
! -------------------------

if (nSym /= 1) then
  call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
end if

! Initializations.
! ----------------

Functional = -9.9d9

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
  write(6,*) SecNam,': RdOne returned ',irc
  write(6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
  call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
end if

if (Debug) then
  write(6,*)
  write(6,*) ' Triangular overlap matrix at start'
  write(6,*) ' ----------------------------------'
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
call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),Name,nBasT,nAtoms,Debug)

! Allocate PA array.
! ------------------
call mma_Allocate(PA,nOrb2LocT,nOrb2LocT,nAtoms,Label='PA')
PA(:,:,:) = 0.0d0

! Localise orbitals.
! ------------------

kOffC = nBasT*nFroT+1
call PipekMezey_Iter(Functional,CMO(kOffC),Work(ipOvlp),Thrs,ThrRot,ThrGrad,PA,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),Name, &
                     nBasT,nOrb2LocT,nAtoms,nMxIter,Maximisation,Converged,Debug,Silent)

! De-allocations.
! ---------------

call mma_deallocate(PA)
call GetMem('nB_Start','Free','Inte',ip_nBas_Start,l_nBas_Start)
call GetMem('nB_per_Atom','Free','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('Ovlp','Free','Real',ipOvlp,lOvlp)

end subroutine PipekMezey
