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

subroutine PipekMezey_Iter(Functional,CMO,Ovlp,Thrs,ThrRot,ThrGrad,PA,nBas_per_Atom,nBas_Start,BName,nBasis,nOrb2Loc,nAtoms, &
                           nMxIter,Maximisation,Converged,Debug,Silent)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nAtoms, nBas_per_Atom(nAtoms), nBas_Start(nAtoms), nBasis, nOrb2Loc, nMxIter
real(kind=wp), intent(out) :: Functional, PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: CMO(nBasis,*)
real(kind=wp), intent(in) :: Ovlp(nBasis,*), Thrs, ThrRot, ThrGrad
character(len=LenIn8), intent(in) :: BName(nBasis)
logical(kind=iwp), intent(in) :: Maximisation, Debug, Silent
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: nIter, i
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2
real(kind=wp), allocatable :: RMat(:,:), PACol(:,:), GradientList(:,:,:), HessianList(:,:),Hdiag_smallList(:,:), &
                            FunctionalList(:)

logical(kind=iwp), parameter :: jacobisweeps = .true.


! Print iteration table header.
! -----------------------------

if (.not. Silent) then
  write(u6,'(//,1X,A,/,1X,A)') '                                                        CPU       Wall', &
                               'nIter       Functional P        Delta     Gradient     (sec)     (sec) %Screen'
end if

! Initialization (iteration 0).
! -----------------------------

if (.not. Silent) call CWTime(C1,W1)
nIter = 0
call mma_Allocate(RMat,nOrb2Loc,nOrb2Loc,Label='RMat')
call mma_Allocate(GradientList,nOrb2Loc,nOrb2Loc,nMxIter,Label='GradientList')  !nMxIter=300, maybe we can make it smaller
call mma_Allocate(FunctionalList,nMxIter,Label='FunctionalList')
FunctionalList(:)=0
call mma_Allocate(HessianList,nOrb2Loc*nOrb2Loc,nMxIter,Label='HessianList')
call mma_Allocate(Hdiag_smallList,nOrb2Loc*(nOrb2Loc+1)/2,nMxIter,Label='Hdiag_smallList')

call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)

FunctionalList(1)=Functional
!write(u6,*) 'In PM_iter: FunctionalList(1) = ', FunctionalList(1)

call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug, GradientList(:,:,1), HessianList(:,1), Hdiag_smallList(:,1))

OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional
if (.not. Silent) then
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,Zero
end if

! Iterations.
! -----------

call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')
Converged = .false.
do while ((nIter < nMxIter) .and. (.not. Converged))
    if (.not. Silent) call CWTime(C1,W1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !choose between optimization methods
    if (jacobisweeps) then
        call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,Maximisation,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,ThrRot,PctSkp,Debug)
    else
        !call newoptimizer(GradientList, FunctionalList)
        call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
    end if

    nIter = nIter+1
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
    FunctionalList(nIter+1)=Functional !first entry is from before first iteration

    if (Debug) then
        write(u6,*) 'nIter = ', nIter
        write(u6,*) 'In PM_iter: FunctionalList(nIter+1) = ', FunctionalList(nIter+1)
    end if

    !calculates nxn gradient matrix for the current iteration and adds it to the List
    call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug,GradientList(:,:,nIter+1), HessianList(:,nIter+1), &
        Hdiag_smallList(:,nIter+1))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !check if converged
    Delta = Functional-OldFunctional
    OldFunctional = Functional
    if (.not. Silent) then
        call CWTime(C2,W2)
        TimC = C2-C1
        TimW = W2-W1
        write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,PctSkp
    end if
    Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)
end do


if (Debug) then
    write(u6,*) 'In PipekMezey_iter'
    write(u6,*) '------------------'
    write(u6,*) 'nIterTot = ', nIter
    do i=1,nIter+1
        write(u6,*) 'for iteration ', i-1
        call RecPrt('GradientList',' ',GradientList(:,:,i), nOrb2Loc, nOrb2Loc)
        call RecPrt('HessianList',' ',HessianList(:,i), nOrb2Loc*nOrb2Loc, 1)
        call RecPrt('Hessian_smallList',' ',Hdiag_smallList(:,i), nOrb2Loc*(nOrb2Loc+1)/2, 1)
    end do
end if

call mma_Deallocate(PACol)
call mma_Deallocate(RMat)
call mma_Deallocate(GradientList)
call mma_Deallocate(FunctionalList)
call mma_Deallocate(HessianList)
call mma_Deallocate(Hdiag_smallList)
! Print convergence message.
! --------------------------

if (.not. Silent) then
    if (.not. Converged) then
        write(u6,'(/,A,I4,A)') 'No convergence after',nIter,' iterations.'
    else
        write(u6,'(/,A,I4,A)') 'Convergence after',nIter,' iterations.'
        write(u6,*)
        write(u6,'(A,I8)') 'Number of localised orbitals  : ',nOrb2loc
        write(u6,'(A,ES20.10)') 'Value of P before localisation: ',FirstFunctional
        write(u6,'(A,ES20.10)') 'Value of P after localisation : ',Functional
    end if
end if

end subroutine PipekMezey_Iter
