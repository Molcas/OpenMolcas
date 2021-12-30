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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine Boys(Functional,CMO,Thrs,ThrRot,ThrGrad,nBas,nOrb2Loc,nFro,nSym,nMxIter,Maximisation,Converged,Debug,Silent)
! Author: T.B. Pedersen
!
! Purpose: Boys localisation of occupied orbitals.

implicit real*8(a-h,o-z)
real*8 CMO(*)
integer nBas(nSym), nOrb2Loc(nSym), nFro(nSym)
logical Maximisation, Converged, Debug, Silent
#include "WrkSpc.fh"
character*4 SecNam
parameter(SecNam='Boys')
parameter(nComp=3) ! 3 components of dipole operator
integer ipLbl(nComp), ipLbl_MO(nComp)
character*8 Label, AlloLbl(nComp), AlloLbl_MO(nComp)

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

! Read AO dipole moment integrals.
! --------------------------------

lLbl = nBasT*nBasT
do iComp=1,nComp
  write(AlloLbl(iComp),'(A,I2)') 'Dipole',iComp
  call GetMem(AlloLbl(iComp),'Allo','Real',ipLbl(iComp),lLbl)
end do

lAux = nBasT*(nBasT+1)/2+4
call GetMem('DipAux','Allo','Real',ipAux,lAux)
Label = 'Mltpl  1'
do iComp=1,nComp
  irc = -1
  iOpt = 2
  iSym = 1
  call RdOne(irc,iOpt,Label,iComp,Work(ipAux),iSym)
  if (irc /= 0) then
    write(6,*) SecNam,': RdOne returned ',irc
    write(6,*) 'Label = ',Label,'   Component = ',iComp
    call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
  end if
  if (Debug) then
    write(6,*)
    write(6,*) ' Triangular dipole matrix at start'
    write(6,*) ' ---------------------------------'
    write(6,*) ' Component: ',iComp
    call TriPrt(' ',' ',Work(ipAux),nBasT)
  end if
  call Tri2Rec(Work(ipAux),Work(ipLbl(iComp)),nBasT,Debug)
end do
call GetMem('DipAux','Free','Real',ipAux,lAux)

! Allocate MO arrays.
! -------------------

lLbl_MO = nOrb2LocT*nOrb2LocT
do iComp=1,nComp
  write(AlloLbl_MO(iComp),'(A,I2)') 'MO dip',iComp
  call GetMem(AlloLbl_MO(iComp),'Allo','Real',ipLbl_MO(iComp),lLbl_MO)
end do

! Localise orbitals.
! ------------------

kOffC = nBasT*nFroT+1
call Boys_Iter(Functional,CMO(kOffC),Thrs,ThrRot,ThrGrad,ipLbl,ipLbl_MO,nBasT,nOrb2LocT,nComp,nMxIter,Maximisation,Converged, &
               Debug,Silent)

! De-allocations.
! ---------------

do iComp=nComp,1,-1
  call GetMem(AlloLbl_MO(iComp),'Free','Real',ipLbl_MO(iComp),lLbl_MO)
end do
do iComp=nComp,1,-1
  call GetMem(AlloLbl(iComp),'Free','Real',ipLbl(iComp),lLbl)
end do

end subroutine Boys
