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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Domain_Localisation(irc)
! Thomas Bondo Pedersen, January 2006.
!
! Purpose: set up orbital domains and pair domains. Find number of
!          strong, weak, distant, and very distant pairs.

use Localisation_globals, only: AnaDomain, BName, CMO, nAtoms, nBas, nFro, nOrb2Loc, nSym, ThrDomain, ThrPairDomain
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
#include "debug.fh"
integer(kind=iwp) :: i, iC, iChange, iCount(0:3), ij, kC, nAtom, nBasT, nnOcc, nOcc
real(kind=wp) :: Fac, ThrPD(3), Tst
integer(kind=iwp), allocatable :: iClass(:), iDomain(:), iPairDomain(:), nBas_per_Atom(:), nBas_Start(:)
real(kind=wp), allocatable :: Coord(:,:), f(:), QD(:), Rmin(:)
character(len=*), parameter :: SecNam = 'Domain_Localisation'

! Set return code.
! ----------------

irc = 0

! Check for symmetry (not allowed).
! ---------------------------------

if (nSym /= 1) then
  irc = -1
  return
end if

! Initializations.
! ----------------

nBasT = nBas(1)
nOcc = nOrb2Loc(1)
nnOcc = nOcc*(nOcc+1)/2
nAtom = nAtoms

! There must be at least 2 atoms and 2 orbitals.
! ----------------------------------------------

if ((nAtom < 2) .or. (nOcc < 2)) then
  irc = -2
  return
end if

! Allocate and get index arrays for indexation of basis functions on
! each atom.
! ------------------------------------------------------------------

call mma_allocate(nBas_per_Atom,nAtom,label='nB_per_Atom')
call mma_allocate(nBas_Start,nAtom,label='nB_Start')
call BasFun_Atom(nBas_per_Atom,nBas_Start,BName,nBasT,nAtom,Debug)

! Define domains.
! ---------------

call mma_allocate(iDomain,(nAtom+1)*nOcc,label='iDomain')
call mma_allocate(QD,nOcc,label='QD')
call mma_allocate(f,nOcc,label='f')

kC = nBasT*nFro(1)+1
call DefineDomain(irc,iDomain,QD,f,CMO(kC),ThrDomain,nBas_per_Atom,nBas_Start,nAtom,nBasT,nOcc)
if (irc /= 0) then
  write(u6,*) SecNam,': ERROR: DefineDomain returned ',irc
  call Error(irc) ! return after deallocations
  return
end if

if (Debug) then
  write(u6,*) SecNam,': checking domain definitions...'
  call CheckDomain(irc,iDomain,nAtom,nOcc)
  if (irc == 0) then
    write(u6,*) '....OK!'
  else
    write(u6,*) '....Ooops. Buggy domain definition!'
    call Error(2) ! return after deallocations
    return
  end if
end if

! Define pair domains.
! Make sure that ThrPairDomain is in ascending order.
! ---------------------------------------------------

ThrPD(:) = ThrPairDomain(:)
call Cho_Order(ThrPD,3,1)
iChange = 0
i = 0
do while ((i < 2) .and. (iChange == 0))
  i = i+1
  Tst = ThrPairDomain(i)-ThrPD(i)
  if (abs(Tst) > 1.0e-15_wp) then
    iChange = 1
  end if
end do

call mma_allocate(iPairDomain,(nAtom+1)*nnOcc,label='iPairDomain')
call mma_allocate(iClass,nnOcc,label='iClass')
call mma_allocate(Rmin,nnOcc,label='Rmin')
call mma_allocate(Coord,3,nAtom,label='Coord')

call Get_dArray('Unique Coordinates',Coord,3*nAtom)
call DefinePairDomain(irc,iPairDomain,iClass,Rmin,iDomain,ThrPD,Coord,nAtom,nOcc,3)
if (irc /= 0) then
  write(u6,*) SecNam,': ERROR: DefinePairDomain returned ',irc
  call Error(irc) ! return after deallocations
  return
end if

if (Debug) then
  write(u6,*) SecNam,': checking pair domain definitions...'
  call CheckDomain(irc,iPairDomain,nAtom,nnOcc)
  if (irc == 0) then
    write(u6,*) '....OK!'
  else
    write(u6,*) '....Ooops. Buggy pair domain definition!'
    call Error(3) ! return after deallocations
    return
  end if
end if

! Print info.
! -----------

call Domain_Histogram(iDomain,nAtom,nOcc,'Histogram of domain sizes')
call Domain_Histogram(iPairDomain,nAtom,nnOcc,'Histogram of pair domain sizes')

call Cho_Head('Pair domain classification','=',80,u6)
iCount(:) = 0
do ij=1,nnOcc
  iC = iClass(ij)
  iCount(iC) = iCount(iC)+1
end do
write(u6,'(/,A)') 'Definition:'
if (iChange /= 0) then
  write(u6,'(A)') 'Notice: the input thresholds were re-ordered to ascending order'
  write(u6,'(A,1P,3(1X,D15.5))') 'Your input order was:',(ThrPairDomain(i),i=1,3)
end if
write(u6,'(A,1P,D15.5)') 'Strong       pairs:                   R <= ',ThrPD(1)
write(u6,'(A,1P,D15.5,A,D15.5)') 'Weak         pairs: ',ThrPD(1),' < R <= ',ThrPD(2)
write(u6,'(A,1P,D15.5,A,D15.5)') 'Distant      pairs: ',ThrPD(2),' < R <= ',ThrPD(3)
write(u6,'(A,1P,D15.5,A)') 'Very distant pairs: ',ThrPD(3),' < R'
write(u6,'(/,A)') 'Classification:'
Fac = 100.0_wp/real(nnOcc,kind=wp)
write(u6,'(A,I9,3X,F7.2,A)') 'Number of strong       pairs: ',iCount(0),Fac*iCount(0),'%'
write(u6,'(A,I9,3X,F7.2,A)') 'Number of weak         pairs: ',iCount(1),Fac*iCount(1),'%'
write(u6,'(A,I9,3X,F7.2,A)') 'Number of distant      pairs: ',iCount(2),Fac*iCount(2),'%'
write(u6,'(A,I9,3X,F7.2,A,/)') 'Number of very distant pairs: ',iCount(3),Fac*iCount(3),'%'

! Analysis of individual domains (if requested).
! ----------------------------------------------

if (AnaDomain) then
  call Analysis_Domain(iDomain,QD,f,Coord,BName,nBas_Start,nAtom,nBasT,nOcc)
end if

! Deallocations.
! --------------

call Error(0)

contains

subroutine Error(code)
  integer(kind=iwp), intent(in) :: code
  if (code /= 0) irc = code
  call mma_deallocate(nBas_per_Atom)
  call mma_deallocate(nBas_Start)
  call mma_deallocate(iDomain)
  call mma_deallocate(QD)
  call mma_deallocate(f)
  if (allocated(iPairDomain)) then
    call mma_deallocate(iPairDomain)
    call mma_deallocate(iClass)
    call mma_deallocate(Rmin)
    call mma_deallocate(Coord)
  end if
end subroutine Error

end subroutine Domain_Localisation
