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

subroutine PolScf(iDist,iDistIm,iDT,iFI,iFP,iFil,iCStart,iTri,VMat,Smat,DiFac,Ract,icnum,Energy,NVarv,iMOC,Haveri,iQ_Atoms, &
                  ip_ExpVal,Poli)

use qmstat_global, only: ChaNuc, CT, DenCorrD, FockM, HHmat, iOcc1, iOrb, iSupM, lExtr, Mp2DensCorr, nPart, nPol, PotNuc, qTot, &
                         Trace_MP2, xyzMyI, xyzMyP, xyzMyQ, xyzQuQ
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iDist, iDistIm, iDT(3), iFI(3), iFP(3), iQ_Atoms, iFil(nTri_Elem(iQ_Atoms),10), iCStart, iTri, icnum, NVarv, &
                     iMOC, ip_ExpVal
real(kind=wp) :: VMat(iTri), Smat(iTri), DiFac, Ract, Energy, Poli(nTri_Elem(iQ_Atoms),10)
logical(kind=iwp) :: Haveri
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iDum, iErr, iGri, iOrba, irr3, iScratch, ixx, ixxi, iyy, iyyi, izz, izzi, j, nFound, nPolCent, nQMCent
real(kind=wp) :: Dummy, Egun, OneEl, PolFac, R2inv, Rinv
logical(kind=iwp) :: JaNej
real(kind=wp), allocatable :: EEigen(:), FFp(:,:), RoMat(:), VpolMat(:)
#include "warnings.h"

! Allocate and initialize the eigenvector matrix with the unit matrix.

iOrba = iOrb(1)
call GetMem('Coeff','Allo','Real',iMOC,iOrba**2)
call dcopy_(iOrba**2,[Zero],0,Work(iMOC),1)
call dcopy_(iOrba,[One],0,Work(iMOC),iOrba+1)

! Define some numbers.

nQMCent = nTri_Elem(iQ_Atoms)
nPolCent = nPart*nPol
Rinv = One/Ract
R2inv = Rinv**2
PolFac = DiFac/Ract**3
Egun = Zero

! Compute some distances needed for the polarization. See polras for some additional information on this.

call Memory_PolPrep('Allo',ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri,nPol,nPart)
call PolPrep(iDist,iDistIM,Work(ixx),Work(iyy),Work(izz),Work(irr3),Work(ixxi),Work(iyyi),Work(izzi),Work(iGri),iCNum,nPolCent)

! Polarization loop commencing.

call mma_allocate(FFp,nPolCent,3,label='FFp')
call mma_allocate(RoMat,max(nTri_Elem(iOrba),iTri),label='RoMat')
call mma_allocate(VpolMat,max(nTri_Elem(iOrba),iTri),label='VpolMat')
call mma_allocate(EEigen,iOrba,label='EEigen')
NVarv = 0
do
  NVarv = NVarv+1
  Energy = Zero
  call PolSolv(iDT,iFI,iFP,Work(ixx),Work(iyy),Work(izz),Work(irr3),Work(ixxi),Work(iyyi),Work(izzi),Work(iGri),FFp,iCNum,r2Inv, &
               DiFac,nPolCent)
  call Densi_MO(Romat,Work(iMOC),1,iOcc1,iOrba,iOrba)
  if (Mp2DensCorr) call DCorrCorr(Romat,DenCorrD,Trace_MP2,iOrba,iOcc1)
  call Polink(Energy,nPolCent,nQMCent,iFil,VpolMat,FFp,PolFac,Poli,iCstart,iTri,iQ_Atoms,qTot,ChaNuc,xyzMyQ,xyzMyI,xyzMyP,Romat, &
              xyzQuQ,CT)

  ! Construct the Fock-matrix from two-electron super-matrix and one-electron matrix, with solvent perturbations added.

  do i=1,iTri
    FockM(i) = Zero
    do j=1,iTri
      FockM(i) = FockM(i)+RoMat(j)*Work(iSupM+iTri*(i-1)+j-1)
    end do
    OneEl = HHmat(i)+Vmat(i)+VpolMat(i)+Smat(i)
    FockM(i) = FockM(i)+OneEl
    ! See Szabo-Ostlund eq. 3.184.
    Energy = Energy+(FockM(i)+OneEl)*RoMat(i)
  end do
  ! Add potential-nuclear energy.
  Energy = Potnuc+Energy*Half

  ! If energy is strange, scream!

  if (Energy > 0) then
    write(u6,*)
    write(u6,*) '  SCF energy is positive. Serious error somewhere.'
    call Quit(_RC_GENERAL_ERROR_)
  end if

  ! Diagonalize the Fock-matrix. Eigenvalues are sorted.

  call GetMem('Scratch','Allo','Real',iScratch,iOrba**2)
  call Diag_Driver('V','A','L',iOrba,FockM,Work(iScratch),iOrba,Dummy,Dummy,iDum,iDum,EEigen,Work(iMOC),iOrba,1,-1,'J',nFound,iErr)
  call GetMem('Scratch','Free','Real',iScratch,iOrba**2)

  ! Check if polarization loop has converged.

  call HaveWeConv(iCNum,iCStart,iQ_Atoms,nPolCent,iDT,FFp,xyzMyI,Egun,Energy,NVarv,JaNej,Haveri)
  if (Haveri .or. JaNej) exit
end do

! Deallocate.

call mma_deallocate(FFp)
call mma_deallocate(RoMat)
call mma_deallocate(EEigen)
call Memory_PolPrep('Free',ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri,nPol,nPart)

! If expectation values are extracted, make a detour.

if (lExtr(6)) call Expectus('SCF  ',HHmat,Vmat,VpolMat,Smat,iMOC,iOrba,.false.,iOcc1,ip_ExpVal)
call mma_deallocate(VpolMat)

! The end is near, hold me!

return

end subroutine PolScf
