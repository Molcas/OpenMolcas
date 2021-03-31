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

use Localisation_globals, only: AnaDomain, ipCMO, BName, nAtoms, nBas, nFro, nOrb2Loc, nSym, ThrDomain, ThrPairDomain
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
#include "WrkSpc.fh"
#include "debug.fh"
integer(kind=iwp) :: i, iC, iChange, iCount(0:3), ij, ip_Coord, ip_f, ip_iClass, ip_iDomain, ip_iPairDomain, ip_nBas_per_Atom, &
                     ip_nBas_Start, ip_QD, ip_Rmin, kC, l_Coord, l_f, l_iClass, l_iDomain, l_iPairDomain, l_nBas_per_Atom, &
                     l_nBas_Start, l_QD, l_Rmin, nAtom, nBasT, nnOcc, nOcc
real(kind=wp) :: Fac, ThrPD(3), Tst
character(len=19), parameter :: SecNam = 'Domain_Localisation'

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

l_nBas_per_Atom = 0
l_nBas_Start = 0
l_iDomain = 0
l_QD = 0
l_f = 0
l_iPairDomain = 0
l_iClass = 0
l_Rmin = 0
l_Coord = 0

! There must be at least 2 atoms and 2 orbitals.
! ----------------------------------------------

if ((nAtom < 2) .or. (nOcc < 2)) then
  irc = -2
  return
end if

! Allocate and get index arrays for indexation of basis functions on
! each atom.
! ------------------------------------------------------------------

l_nBas_per_Atom = nAtom
l_nBas_Start = nAtom
call GetMem('nB_per_Atom','Allo','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('nB_Start','Allo','Inte',ip_nBas_Start,l_nBas_Start)
call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),BName,nBasT,nAtom,Debug)

! Define domains.
! ---------------

l_iDomain = (nAtom+1)*nOcc
l_QD = nOcc
l_f = nOcc
call GetMem('iDomain','Allo','Inte',ip_iDomain,l_iDomain)
call GetMem('QD','Allo','Real',ip_QD,l_QD)
call GetMem('f','Allo','Real',ip_f,l_f)

kC = ipCMO+nBasT*nFro(1)
call DefineDomain(irc,iWork(ip_iDomain),Work(ip_QD),Work(ip_f),Work(kC),ThrDomain,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start), &
                  nAtom,nBasT,nOcc)
if (irc /= 0) then
  write(u6,*) SecNam,': ERROR: DefineDomain returned ',irc
  call Error(irc) ! return after deallocations
  return
end if

if (Debug) then
  write(u6,*) SecNam,': checking domain definitions...'
  call CheckDomain(irc,iWork(ip_iDomain),nAtom,nOcc)
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

call dCopy_(3,ThrPairDomain,1,ThrPD,1)
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

l_iPairDomain = (nAtom+1)*nnOcc
l_iClass = nnOcc
l_Rmin = nnOcc
l_Coord = 3*nAtom
call GetMem('iPairDomain','Allo','Inte',ip_iPairDomain,l_iPairDomain)
call GetMem('iClass','Allo','Inte',ip_iClass,l_iClass)
call GetMem('Rmin','Allo','Real',ip_Rmin,l_Rmin)
call GetMem('Coord','Allo','Real',ip_Coord,l_Coord)

call Get_dArray('Unique Coordinates',Work(ip_Coord),l_Coord)
call DefinePairDomain(irc,iWork(ip_iPairDomain),iWork(ip_iClass),Work(ip_Rmin),iWork(ip_iDomain),ThrPD,Work(ip_Coord),nAtom,nOcc,3)
if (irc /= 0) then
  write(u6,*) SecNam,': ERROR: DefinePairDomain returned ',irc
  call Error(irc) ! return after deallocations
  return
end if

if (Debug) then
  write(u6,*) SecNam,': checking pair domain definitions...'
  call CheckDomain(irc,iWork(ip_iPairDomain),nAtom,nnOcc)
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

call Domain_Histogram(iWork(ip_iDomain),nAtom,nOcc,'Histogram of domain sizes')
call Domain_Histogram(iWork(ip_iPairDomain),nAtom,nnOcc,'Histogram of pair domain sizes')

call Cho_Head('Pair domain classification','=',80,u6)
do i=0,3
  iCount(i) = 0
end do
do ij=0,nnOcc-1
  iC = iWork(ip_iClass+ij)
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
  call Analysis_Domain(iWork(ip_iDomain),Work(ip_QD),Work(ip_f),Work(ip_Coord),BName,iWork(ip_nBas_Start),nAtom,nBasT,nOcc)
end if

! Deallocations.
! --------------

call Error(0)

contains

subroutine Error(code)
  integer(kind=iwp), intent(in) :: code
  if (code /= 0) irc = code
  if (l_Coord > 0) then
    call GetMem('Coord','Free','Real',ip_Coord,l_Coord)
  end if
  if (l_Rmin > 0) then
    call GetMem('Rmin','Free','Real',ip_Rmin,l_Rmin)
  end if
  if (l_iClass > 0) then
    call GetMem('iClass','Free','Inte',ip_iClass,l_iClass)
  end if
  if (l_iPairDomain > 0) then
    call GetMem('iPairDomain','Free','Inte',ip_iPairDomain,l_iPairDomain)
  end if
  if (l_f > 0) then
    call GetMem('f','Free','Real',ip_f,l_f)
  end if
  if (l_QD > 0) then
    call GetMem('QD','Free','Real',ip_QD,l_QD)
  end if
  if (l_iDomain > 0) then
    call GetMem('iDomain','Free','Inte',ip_iDomain,l_iDomain)
  end if
  if (l_nBas_Start > 0) then
    call GetMem('nB_Start','Free','Inte',ip_nBas_Start,l_nBas_Start)
  end if
  if (l_nBas_per_Atom > 0) then
    call GetMem('nB_per_Atom','Free','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
  end if
end subroutine Error

end subroutine Domain_Localisation
