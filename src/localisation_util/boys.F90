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

subroutine Boys(Functional,CMO,nBas,nOrb2Loc,nFro,nSym,nMxIter,Converged,Silent)
! Author: T.B. Pedersen
!
! Purpose: Boys localisation of occupied orbitals.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug

implicit none
real(kind=wp), intent(out) :: Functional
real(kind=wp), intent(inout) :: CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb2Loc(nSym), nFro(nSym), nMxIter
logical(kind=iwp), intent(in) :: Silent
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: iCmp, iComp, iOpt, irc, iSym, kOffC, lAux, nBasT, nFroT, nOrb2LocT
character(len=8) :: Label
real(kind=wp), allocatable :: Aux(:), Lbl(:,:,:), Lbl_AO(:,:,:), SMat(:)
real(kind=wp) :: CoM(3)
integer(kind=iwp), parameter :: nComp = 3 ! 3 components of dipole operator
character(len=*), parameter :: SecNam = 'Boys'

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

! Read AO dipole moment integrals.
! --------------------------------

call mma_allocate(Lbl_AO,nBasT,nBasT,nComp,label='Dipole')

lAux = nBasT*(nBasT+1)/2+4

call mma_allocate(SMat,lAux,label='SMat')
Label = 'Mltpl  0'
iCmp = 1
irc = -1
iOpt = 0
iSym = 1
call RdOne(irc,iOpt,Label,iCmp,SMat,iSym)
Call Get_dArray('Center of Mass',CoM,3)

call mma_allocate(Aux,lAux,label='DipAux')
Label = 'Mltpl  1'
do iComp=1,nComp
  iCmp = iComp
  irc = -1
  iOpt = 0
  iSym = 1
  call RdOne(irc,iOpt,Label,iCmp,Aux,iSym)
  if (irc /= 0) then
    write(u6,*) SecNam,': RdOne returned ',irc
    write(u6,*) 'Label = ',Label,'   Component = ',iComp
    call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
  end if
  if (Debug) then
    write(u6,*)
    write(u6,*) ' Triangular dipole matrix at start'
    write(u6,*) ' ---------------------------------'
    write(u6,*) ' Component: ',iComp
    call TriPrt(' ',' ',Aux,nBasT)
  end if
  Aux(:)=Aux(:) + CoM(iComp)*SMat(:)
  call Tri2Rec(Aux,Lbl_AO(:,:,iComp),nBasT)
end do
call mma_deallocate(SMat)
call mma_deallocate(Aux)

! Allocate MO arrays.
! -------------------

call mma_allocate(Lbl,nOrb2LocT,nOrb2LocT,nComp,label='MO_dip')

! Localise orbitals.
! ------------------

kOffC = 1+nBasT*nFroT
call Boys_Iter(Functional,CMO(kOffC),Lbl_AO,Lbl,nBasT,nOrb2LocT,nComp,nMxIter,Converged, &
               Silent)

! De-allocations.
! ---------------

call mma_deallocate(Lbl_AO)
call mma_deallocate(Lbl)

end subroutine Boys
