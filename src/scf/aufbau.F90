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

use InfSCF, only: nFro, nOcc, nOrb, nSym, rTemp, TEEE
use SCF_Arrays, only: EOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: nAuf(2), nOccup, iOK, nD
real(kind=wp) :: Occup(nOccup,nD)
#include "Molcas.fh"
integer(kind=iwp) :: iD, iOrb, iOrBas, ipOcc, iSym, jOrBas, jSym, mD, mOrb_AS(2), nElec, nEOrb, nOrb_AS(2), nOrBas, Tmp
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

iOrbAS = 1
do iSym=1,nSym
  do iOrb=1,nOrb(iSym)-nFro(iSym)
    do iD=1,nD
      Irp(iOrbAS,iD) = iSym
      Map(iOrbAS,iD) = iOrbAS
    end do
    iOrbAS = iOrbAS+1
  end do
end do
nOrbAS = iOrbAS-1

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
call dcopy_(nOccup*nD,[Zero],0,Occup,1)

if (Teee) then

  UHF_occ = Three-real(nD,kind=wp)
  mD = 2/nD
  do iD=1,nD
    eferm = FermiPop(EOrb(1,iD),Occup(1,iD),nOrbAS,RTemp,nAuf(iD)*mD,UHF_occ)
#   ifdef _DEBUGPRINT_
    write(u6,'(A,G20.10)') '         E(Fermi)=',eferm
#   endif
  end do

  iOrbAS = 0
  do iSym=1,nSym
    do iD=1,nD
      nOrb_AS(iD) = 0
      mOrb_AS(iD) = 0
      sum_el(iD) = Zero
    end do

    jOrbAS = iOrbAS
    unlikelyOcc = 0.19_wp
    do iOrb=1,nOrb(iSym)-nFro(iSym)
      iOrbAS = iOrbAS+1
      do iD=1,nD
        if (Occup(iOrbAS,iD) >= UHF_occ-unlikelyOcc) nOrb_AS(iD) = nOrb_AS(iD)+1
        if (Occup(iOrbAS,iD) < unlikelyOcc) mOrb_AS(iD) = mOrb_AS(iD)+1
        sum_el(iD) = sum_el(iD)+Occup(iOrbAS,iD)
      end do
    end do
    Fact = nD*half
    Fact2 = 0.99_wp+real(2-nD,kind=wp)
    do iD=1,nD
      nOccAuf(iSym,kOccAuf,iD) = nOrb_AS(iD)
      nOcc(iSym,iD) = int(Fact*(sum_el(iD)+Fact2/nSym))
    end do
  end do
  kOccAuf = 3-kOccAuf

else

  Fact = Two/real(nD,kind=wp)
  do iD=1,nD
    do iOrbAS=1,nAuf(iD)
      iSym = Irp(Map(iOrbAS,iD),iD)
      nOcc(iSym,iD) = nOcc(iSym,iD)+1
      ipOcc = 0
      do jSym=1,iSym-1
        ipOcc = ipOcc+nOrb(jSym)
      end do
      Occup(ipOcc+nOcc(iSym,iD),iD) = Fact
    end do
  end do

end if

iOK = 1
do iD=1,nD
  nElec = 0
  do iSym=1,nSym
    if (nOccAuf(iSym,1,iD) /= nOccAuf(iSym,2,iD)) iOK = 0
    nElec = nElec+nOccAuf(iSym,1,iD)
  end do
  if (nElec /= nAuf(iD)) iOK = 0
end do
if (iOK == 1) then
  do iSym=1,nSym
    do iD=1,nD
      nOcc(iSym,iD) = nOccAuf(iSym,1,iD)
    end do
  end do
end if

! Write new occupation on the RUNFILE

call Put_iArray('nIsh',nOcc(1,1),nSym)
if (nD == 2) call Put_iArray('nIsh beta',nOcc(1,2),nSym)

call mma_deallocate(Irp)
call mma_deallocate(Map)

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_real(eferm)
#endif

end subroutine Aufbau
