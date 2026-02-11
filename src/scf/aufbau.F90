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
! Copyright (C) 1995, Martin Schuetz                                   *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Aufbau(nAuf,Occup,nOccup,iOK,nD)
!***********************************************************************
!                                                                      *
!     purpose: sets the orbital occupation numbers in the different    *
!              irreps according to the Aufbau scheme...                *
!                                                                      *
!     method:  sets up a map vector, which is sorted with respect to   *
!              the orbital energies, neglecting the boundaries of the  *
!              different irrep blocks. The lowest orbitals then are    *
!              occupied...                                             *
!     input:                                                           *
!       nAuf          : # (doubly) occupied orbitals                   *
!                                                                      *
!     output:                                                          *
!       Occup(nOccup) : orbital occupation numbers                     *
!***********************************************************************

use InfSCF, only: EOrb, nFro, nOcc, nOrb, nSym, rTemp, TEEE
use Molcas, only: MxSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nAuf(2), nOccup, nD
real(kind=wp), intent(out) :: Occup(nOccup,nD)
integer(kind=iwp), intent(out) :: iOK
integer(kind=iwp) :: iD, iOrb, iOrBas, ipOcc, iSym, jOrBas, mD, mOrb_AS(2), nElec, nEOrb, nOrb_AS(2), nOrBas, Tmp
! These occupation number vectors are used to determine if we have convergence.
integer(kind=iwp) :: kOccAuf = 1, nOccAuf(MxSym,2,2) = -1
real(kind=wp) :: EFerm, Fact, Fact2, Sum_el(2), UHF_Occ, UnlikelyOcc
integer(kind=iwp), allocatable :: Map(:,:), Irp(:,:)
real(kind=wp), external :: Fermipop

nEOrb = size(EOrb,1)

call mma_allocate(Map,nEOrb,nD,Label='Map')
call mma_allocate(Irp,nEOrb,nD,Label='Irp')
!----------------------------------------------------------------------*
! Initialize convergence detection                                     *
!----------------------------------------------------------------------*
if (kOccAuf == -1) then
  nOccAuf(:,:,:) = -1
  kOccAuf = 1
end if

! Set up map...

iOrbAS = 0
do iSym=1,nSym
  do iOrb=1,nOrb(iSym)-nFro(iSym)
    iOrbAS = iOrbAS+1
    Irp(iOrbAS,:) = iSym
    Map(iOrbAS,:) = iOrbAS
  end do
end do
nOrbAS = iOrbAS

! Now sort map with respect to orbital energies (bubblesort)

do iOrbAS=1,nOrbAS-1
  do jOrbAS=nOrbAS-1,iOrbAS,-1
    do iD=1,nD
      if (EOrb(Map(jOrbAS,iD),iD) > EOrb(Map(1+jOrbAS,iD),iD)) then
        Tmp = Map(jOrbAS,iD)
        Map(jOrbAS,iD) = Map(1+jOrbAS,iD)
        Map(1+jOrbAS,iD) = Tmp
      end if
    end do
  end do
end do

! and fill up the orbitals...

nOcc(:,:) = 0
Occup(:,:) = Zero

if (Teee) then

  UHF_occ = Three-real(nD,kind=wp)
  mD = 2/nD
  do iD=1,nD
    eferm = FermiPop(EOrb(:,iD),Occup(:,iD),nOrbAS,RTemp,nAuf(iD)*mD,UHF_occ)
#   ifdef _DEBUGPRINT_
    write(u6,'(A,G20.10)') '         E(Fermi)=',eferm
#   else
#   include "macros.fh"
    unused_var(eferm)
#   endif
  end do

  iOrbAS = 0
  do iSym=1,nSym
    nOrb_AS(1:nD) = 0
    mOrb_AS(1:nD) = 0
    sum_el(:) = Zero

    jOrbAS = iOrbAS
    unlikelyOcc = 0.19_wp
    do iOrb=1,nOrb(iSym)-nFro(iSym)
      iOrbAS = iOrbAS+1
      do iD=1,nD
        if (Occup(iOrbAS,iD) >= UHF_occ-unlikelyOcc) nOrb_AS(iD) = nOrb_AS(iD)+1
        if (Occup(iOrbAS,iD) < unlikelyOcc) mOrb_AS(iD) = mOrb_AS(iD)+1
      end do
      sum_el(1:nD) = sum_el(1:nD)+Occup(iOrbAS,1:nD)
    end do
    Fact = nD*Half
    Fact2 = 0.99_wp+real(2-nD,kind=wp)
    nOccAuf(iSym,kOccAuf,1:nD) = nOrb_AS(1:nD)
    nOcc(iSym,1:nD) = int(Fact*(sum_el(1:nD)+Fact2/nSym))
  end do
  kOccAuf = 3-kOccAuf

else

  Fact = Two/real(nD,kind=wp)
  do iD=1,nD
    do iOrbAS=1,nAuf(iD)
      iSym = Irp(Map(iOrbAS,iD),iD)
      nOcc(iSym,iD) = nOcc(iSym,iD)+1
      ipOcc = sum(nOrb(1:iSym-1))
      Occup(ipOcc+nOcc(iSym,iD),iD) = Fact
    end do
  end do

end if

iOK = 1
do iD=1,nD
  if (any(nOccAuf(1:nSym,1,iD) /= nOccAuf(1:nSym,2,iD))) iOK = 0
  nElec = sum(nOccAuf(1:nSym,1,iD))
  if (nElec /= nAuf(iD)) iOK = 0
end do
if (iOK == 1) nOcc(1:nSym,1:nD) = nOccAuf(1:nSym,1,1:nD)

! Write new occupation on the RUNFILE

call Put_iArray('nIsh',nOcc(1,1),nSym)
if (nD == 2) call Put_iArray('nIsh beta',nOcc(1,2),nSym)

call mma_deallocate(Irp)
call mma_deallocate(Map)

return

end subroutine Aufbau
