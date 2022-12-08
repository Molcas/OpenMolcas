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

subroutine PolRas(Dist,DistIm,DT,FI,FP,Fil,iCStart,iTriState,VMat,Smat,DiFac,Ract,icnum,Energy,NVarv,STC,Haveri,iQ_Atoms,ExpVal, &
                  Poli)

use qmstat_global, only: ChaNuc, CT, dCIRef, HmatSOld, HmatState, iCIInd, iExtr_Eig, lCiSelect, lExtr, nCent, nCIRef, nEqState, &
                         nPart, nPol, nState, qTot, xyzMyI, xyzMyP, xyzMyQ, xyzQuQ
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms, iCStart, iTriState, icnum
real(kind=wp), intent(in) :: Dist(nCent,nCent,nTri_Elem(nPart-icnum-1)), DistIm(nCent,nPart-icnum,nCent,nPart-icnum), &
                             FP(3,nPol*nPart), Fil(nPol*nPart,3,nTri_Elem(iQ_Atoms),10), VMat(iTriState), Smat(iTriState), DiFac, &
                             Ract
real(kind=wp), intent(inout) :: DT(3,nPol*nPart), Energy
real(kind=wp), intent(out) :: FI(3,nPol*nPart), STC(nState,nState), ExpVal(4,nState), Poli(nTri_Elem(iQ_Atoms),10)
integer(kind=iwp), intent(out) :: NVarv
logical(kind=iwp), intent(out) :: Haveri
integer(kind=iwp) :: iDum, iErr, nFound, nPolCent, nQMCent
real(kind=wp) :: Dummy, Egun, PolFac, R2inv, Rinv
logical(kind=iwp) :: JaNej
real(kind=wp), allocatable :: EEigen(:), FFp(:,:), Gri(:,:), RoMatSt(:), rr3(:,:), Scratch(:,:), VpolMat(:), xx(:,:), xxi(:,:), &
                              yy(:,:), yyi(:,:), zz(:,:), zzi(:,:)

! Allocate and initialize the eigenvector matrix with the unit matrix.

call unitmat(STC,nState)

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
call mma_allocate(RoMatSt,nTri_Elem(nState),label='RoMatSt')
call mma_allocate(VpolMat,nTri_Elem(nState),label='VpolMat')
call mma_allocate(EEigen,nState,label='EEigen')
NVarv = 0
do
  NVarv = NVarv+1
  Energy = Zero
  call PolSolv(DT,FI,FP,xx,yy,zz,rr3,xxi,yyi,zzi,Gri,FFp,iCNum,r2Inv,DiFac,nPolCent)
  call DensiSt(RomatSt,STC,nEqState,nState,nState)
  call Polins(Energy,nPolCent,nQMCent,Fil,VpolMat,FFp,PolFac,poli,xyzMyQ,xyzMyI,xyzMyP,iCstart,iQ_Atoms,qTot,ChaNuc,RoMatSt, &
              xyzQuQ,CT)

  ! Assemble the Hamiltonian matrix.

  HmatState(:) = HmatSOld+Vmat+VpolMat+SMat
  Energy = Half*Energy

  ! Diagonalize the bastard. Eigenvalues are sorted and the relevant eigenvalue is added to total energy.

  call mma_allocate(Scratch,nState,nState,label='Scratch')
  call Diag_Driver('V','A','L',nState,HMatState,Scratch,nState,Dummy,Dummy,iDum,iDum,EEigen,STC,nState,1,-1,'J',nFound,iErr)
  if (lCiSelect) call CiSelector(nEqState,nState,STC,nCIRef,iCIInd,dCIRef)
  Energy = Energy+EEigen(nEqState)
  call mma_deallocate(Scratch)

  ! Check if polarization loop has converged.

  call HaveWeConv(iCNum,iCStart,iQ_Atoms,nPolCent,DT,FFp,xyzMyI,Egun,Energy,NVarv,JaNej,Haveri)
  if (Haveri .or. JaNej) exit
end do

! Deallocate.

call mma_deallocate(FFp)
call mma_deallocate(RoMatSt)
call mma_deallocate(VpolMat)
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

if (lExtr(6)) call Expectus('RASSI',HmatSOld,Vmat,VpolMat,Smat,STC,nState,lExtr(4),iExtr_Eig,ExpVal)

! Is it dead? It's terminated!

return

end subroutine PolRas
