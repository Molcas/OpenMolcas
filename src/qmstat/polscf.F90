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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm1.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.h"
dimension Poli(MxQCen,10), VMat(MxOT), FFp(nPol*nPart,3)
dimension VpolMat(MxOT), Smat(MxOT), RoMat(MxOT)
dimension EEigen(MxOrb)
dimension iDT(3), iFI(3), iFP(3), iFil(MxQCen,10)
logical JaNej, Haveri

! Allocate and initialize the eigenvector matrix with the unit matrix.

iOrba = iOrb(1)
call GetMem('Coeff','Allo','Real',iMOC,iOrba**2)
call dcopy_(iOrba**2,[ZERO],iZERO,Work(iMOC),iONE)
call dcopy_(iOrba,[ONE],iZERO,Work(iMOC),iOrba+1)

! Define some numbers.

nQMCent = (iQ_Atoms*(iQ_Atoms+1))/2
nPolCent = nPart*nPol
Rinv = 1.0d0/Ract
R2inv = Rinv**2
PolFac = DiFac/Ract**3
Egun = 0.0d0

! Compute some distances needed for the polarization. See polras for some additional information on this.

call Memory_PolPrep('Allo',ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri,nPol,nPart)
call PolPrep(iDist,iDistIM,Work(ixx),Work(iyy),Work(izz),Work(irr3),Work(ixxi),Work(iyyi),Work(izzi),Work(iGri),iCNum,nPolCent)

! Polarization loop commencing.

NVarv = 0
do
  NVarv = NVarv+1
  Energy = 0.0d0
  call PolSolv(iDT,iFI,iFP,Work(ixx),Work(iyy),Work(izz),Work(irr3),Work(ixxi),Work(iyyi),Work(izzi),Work(iGri),FFp,iCNum,r2Inv, &
               DiFac,nPolCent)
  call Densi_MO(Romat,Work(iMOC),1,iOcc1,iOrb(1),iOrb(1))
  if (Mp2DensCorr) call DCorrCorr(Romat,DenCorrD,Trace_MP2,iOrb(1),iOcc1)
  call Polink(Energy,iCall,nPolCent,nQMCent,iFil,VpolMat,FFp,PolFac,Poli,iCstart,iTri,iQ_Atoms,qTot,ChaNuc,xyzMyQ,xyzMyI,xyzMyP, &
              Romat,xyzQuQ,CT)

  ! Construct the Fock-matrix from two-electron super-matrix and one-electron matrix, with solvent perturbations added.

  do i=1,iTri
    FockM(i) = 0.0d0
    do j=1,iTri
      FockM(i) = FockM(i)+RoMat(j)*Work(iSupM+iTri*(i-1)+j-1)
    end do
    OneEl = HHmat(i)+Vmat(i)+VpolMat(i)+Smat(i)
    FockM(i) = FockM(i)+OneEl
    ! See Szabo-Ostlund eq. 3.184.
    Energy = Energy+(FockM(i)+OneEl)*RoMat(i)
  end do
  ! Add potential-nuclear energy.
  Energy = Potnuc+Energy*0.5

  ! If energy is strange, scream!

  if (Energy > 0) then
    write(6,*)
    write(6,*) '  SCF energy is positive. Serious error somewhere.'
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

call Memory_PolPrep('Free',ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri,nPol,nPart)

! If expectation values are extracted, make a detour.

if (lExtr(6)) call Expectus('SCF  ',HHmat,Vmat,VpolMat,Smat,MxOT,iMOC,iOrba,.false.,iOcc1,ip_ExpVal)

! The end is near, hold me!

return

end subroutine PolScf
