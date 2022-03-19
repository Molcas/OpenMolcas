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

subroutine PolRas(iDist,iDistIm,iDT,iFI,iFP,iFil,iCStart,iTriState,VMat,Smat,DiFac,Ract,icnum,Energy,NVarv,iSTC,Haveri,iQ_Atoms, &
                  ip_ExpVal,Poli)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iDist, iDistIm, iDT(3), iFI(3), iFP(3), iQ_Atoms, iFil(nTri_Elem(iQ_Atoms),10), iCStart, iTriState, icnum, &
                     NVarv, iSTC, ip_ExpVal
real(kind=wp) :: VMat(iTriState), Smat(iTriState), DiFac, Ract, Energy, Poli(nTri_Elem(iQ_Atoms),10)
logical(kind=iwp) :: Haveri
integer(kind=iwp) :: i, iDum, iErr, iGri, irr3, iScratch, ixx, ixxi, iyy, iyyi, izz, izzi, nFound, nPolCent, nQMCent
real(kind=wp) :: Dummy, Egun, PolFac, R2inv, Rinv
logical(kind=iwp) :: JaNej
real(kind=wp), allocatable :: EEigen(:), FFp(:,:), RoMatSt(:), VpolMat(:)

! Allocate and initialize the eigenvector matrix with the unit matrix.

call GetMem('Coeff','Allo','Real',iSTC,nState**2)
call dcopy_(nState**2,[Zero],0,Work(iSTC),1)
call dcopy_(nState,[One],0,Work(iSTC),nState+1)

! Define some numbers.

nQMCent = nTri_Elem(iQ_Atoms)
nPolCent = nPart*nPol
Rinv = One/Ract
R2inv = Rinv**2
PolFac = DiFac/Ract**3
Egun = Zero

! Compute distances and vectors needed for solving the polarization
! equations. Observe the use of a memory allocator before the call
! to polprep. This is to conserve memory but without having to
! rewrite the entire polprep routine (written originally with
! static allocations). A future project is to rewrite polprep.

call Memory_PolPrep('Allo',ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri,nPol,nPart)
call PolPrep(iDist,iDistIM,Work(ixx),Work(iyy),Work(izz),Work(irr3),Work(ixxi),Work(iyyi),Work(izzi),Work(iGri),iCNum,nPolCent)

! Polarization loop commencing.

call mma_allocate(FFp,nPolCent,3,label='FFp')
call mma_allocate(RoMatSt,nTri_Elem(nState),label='RoMatSt')
call mma_allocate(VpolMat,nTri_Elem(nState),label='VpolMat')
call mma_allocate(EEigen,nState,label='EEigen')
NVarv = 0
do
  NVarv = NVarv+1
  Energy = Zero
  call PolSolv(iDT,iFI,iFP,Work(ixx),Work(iyy),Work(izz),Work(irr3),Work(ixxi),Work(iyyi),Work(izzi),Work(iGri),FFp,iCNum,r2Inv, &
               DiFac,nPolCent)
  call DensiSt(RomatSt,Work(iSTC),nEqState,nState,nState)
  call Polins(Energy,nPolCent,nQMCent,iFil,VpolMat,FFp,PolFac,poli,xyzMyQ,xyzMyI,xyzMyP,iCstart,iQ_Atoms,qTot,ChaNuc,RoMatSt, &
              xyzQuQ,CT)

  ! Assemble the Hamiltonian matrix.

  do i=1,iTriState
    HmatState(i) = HmatSOld(i)+Vmat(i)+VpolMat(i)+SMat(i)
  end do
  Energy = Half*Energy

  ! Diagonalize the bastard. Eigenvalues are sorted and the relevant eigenvalue is added to total energy.

  call GetMem('Scratch','Allo','Real',iScratch,nState**2)
  call Diag_Driver('V','A','L',nState,HMatState,Work(iScratch),nState,Dummy,Dummy,iDum,iDum,EEigen,Work(iSTC),nState,1,-1,'J', &
                   nFound,iErr)
  if (lCiSelect) call CiSelector(nEqState,nState,iSTC,nCIRef,iCIInd,dCIRef)
  Energy = Energy+EEigen(nEqState)
  call GetMem('Scratch','Free','Real',iScratch,nState**2)

  ! Check if polarization loop has converged.

  call HaveWeConv(iCNum,iCStart,iQ_Atoms,nPolCent,iDT,FFp,xyzMyI,Egun,Energy,NVarv,JaNej,Haveri)
  if (Haveri .or. JaNej) exit
end do

! Deallocate.

call mma_deallocate(FFp)
call mma_deallocate(RoMatSt)
call mma_deallocate(VpolMat)
call mma_deallocate(EEigen)
call Memory_PolPrep('Free',ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri,nPol,nPart)

! If expectation values are extracted, make a detour.

if (lExtr(6)) call Expectus('RASSI',HmatSOld,Vmat,VpolMat,Smat,iSTC,nState,lExtr(4),iExtr_Eig,ip_ExpVal)

! Is it dead? It's terminated!

return

end subroutine PolRas
