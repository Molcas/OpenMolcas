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
integer(kind=iwp) :: nIter
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2
real(kind=wp), allocatable :: RMat(:,:), PACol(:,:)

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
call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug)
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
  call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,Maximisation,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,ThrRot,PctSkp,Debug)
  call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
  call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug)
  nIter = nIter+1
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
call mma_Deallocate(PACol)
call mma_Deallocate(RMat)

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
