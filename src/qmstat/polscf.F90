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

subroutine PolScf(Dist,DistIm,DT,FI,FP,Fil,iCStart,iTri,VMat,Smat,DiFac,Ract,icnum,Energy,NVarv,MOC,Haveri,iQ_Atoms,ExpVal,Poli)

use qmstat_global, only: ChaNuc, CT, DenCorrD, FockM, HHmat, iOcc1, iOrb, lExtr, Mp2DensCorr, nCent, nPart, nPol, PotNuc, qTot, &
                         SupM, Trace_MP2, xyzMyI, xyzMyP, xyzMyQ, xyzQuQ
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms, iCStart, iTri, icnum
real(kind=wp), intent(in) :: Dist(nCent,nCent,nTri_Elem(nPart-icnum-1)), DistIm(nCent,nPart-icnum,nCent,nPart-icnum), &
                             FP(3,nPol*nPart), Fil(nPol*nPart,3,nTri_Elem(iQ_Atoms),10), VMat(iTri), Smat(iTri), DiFac, Ract
real(kind=wp), intent(inout) :: DT(3,nPol*nPart), Energy
real(kind=wp), intent(out) :: FI(3,nPol*nPart), MOC(iOrb(1),iOrb(1)), ExpVal(4,1), Poli(nTri_Elem(iQ_Atoms),10)
integer(kind=iwp), intent(out) :: NVarv
logical(kind=iwp), intent(out) :: Haveri
integer(kind=iwp) :: i, iDum, iErr, iOrba, j, nFound, nPolCent, nQMCent
real(kind=wp) :: Dummy, Egun, OneEl, PolFac, R2inv, Rinv
logical(kind=iwp) :: JaNej
real(kind=wp), allocatable :: EEigen(:), FFp(:,:), Gri(:,:), RoMat(:), rr3(:,:), Scratch(:,:), VpolMat(:), xx(:,:), xxi(:,:), &
                              yy(:,:), yyi(:,:), zz(:,:), zzi(:,:)
#include "warnings.h"

! Allocate and initialize the eigenvector matrix with the unit matrix.

iOrba = iOrb(1)
call unitmat(MOC,iOrba)

! Define some numbers.

nQMCent = nTri_Elem(iQ_Atoms)
nPolCent = nPart*nPol
Rinv = One/Ract
R2inv = Rinv**2
PolFac = DiFac/Ract**3
Egun = Zero

! Compute some distances needed for the polarization. See polras for some additional information on this.

call mma_allocate(xx,nPolCent,nPolCent,label='xx')
call mma_allocate(yy,nPolCent,nPolCent,label='yy')
call mma_allocate(zz,nPolCent,nPolCent,label='zz')
call mma_allocate(xxi,nPolCent,nPolCent,label='ixx')
call mma_allocate(yyi,nPolCent,nPolCent,label='iyy')
call mma_allocate(zzi,nPolCent,nPolCent,label='izz')
call mma_allocate(rr3,nPolCent,nPolCent,label='irr3')
call mma_allocate(Gri,nPolCent,nPolCent,label='iGri')
xx(:,:) = Zero
yy(:,:) = Zero
zz(:,:) = Zero
xxi(:,:) = Zero
yyi(:,:) = Zero
zzi(:,:) = Zero
rr3(:,:) = Zero
Gri(:,:) = Zero
call PolPrep(Dist,DistIm,xx,yy,zz,rr3,xxi,yyi,zzi,Gri,iCNum,nPolCent)

! Polarization loop commencing.

call mma_allocate(FFp,nPolCent,3,label='FFp')
call mma_allocate(RoMat,max(nTri_Elem(iOrba),iTri),label='RoMat')
call mma_allocate(VpolMat,max(nTri_Elem(iOrba),iTri),label='VpolMat')
call mma_allocate(EEigen,iOrba,label='EEigen')
NVarv = 0
do
  NVarv = NVarv+1
  Energy = Zero
  call PolSolv(DT,FI,FP,xx,yy,zz,rr3,xxi,yyi,zzi,Gri,FFp,iCNum,r2Inv,DiFac,nPolCent)
  call Densi_MO(Romat,MOC,1,iOcc1,iOrba,iOrba)
  if (Mp2DensCorr) call DCorrCorr(Romat,DenCorrD,Trace_MP2,iOrba,iOcc1)
  call Polink(Energy,nPolCent,nQMCent,Fil,VpolMat,FFp,PolFac,Poli,iCstart,iTri,iQ_Atoms,qTot,ChaNuc,xyzMyQ,xyzMyI,xyzMyP,Romat, &
              xyzQuQ,CT)

  ! Construct the Fock-matrix from two-electron super-matrix and one-electron matrix, with solvent perturbations added.

  FockM(:) = Zero
  do i=1,iTri
    do j=1,iTri
      FockM(i) = FockM(i)+RoMat(j)*SupM(j,i)
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

  call mma_allocate(Scratch,iOrba,iOrba,label='Scratch')
  call Diag_Driver('V','A','L',iOrba,FockM,Scratch,iOrba,Dummy,Dummy,iDum,iDum,EEigen,MOC,iOrba,1,-1,'J',nFound,iErr)
  call mma_deallocate(Scratch)

  ! Check if polarization loop has converged.

  call HaveWeConv(iCNum,iCStart,iQ_Atoms,nPolCent,DT,FFp,xyzMyI,Egun,Energy,NVarv,JaNej,Haveri)
  if (Haveri .or. JaNej) exit
end do

! Deallocate.

call mma_deallocate(FFp)
call mma_deallocate(RoMat)
call mma_deallocate(EEigen)
call mma_deallocate(xx)
call mma_deallocate(yy)
call mma_deallocate(zz)
call mma_deallocate(xxi)
call mma_deallocate(yyi)
call mma_deallocate(zzi)
call mma_deallocate(rr3)
call mma_deallocate(Gri)

! If expectation values are extracted, make a detour.

if (lExtr(6)) call Expectus('SCF  ',HHmat,Vmat,VpolMat,Smat,MOC,iOrba,.false.,iOcc1,ExpVal)
call mma_deallocate(VpolMat)

! The end is near, hold me!

return

end subroutine PolScf
