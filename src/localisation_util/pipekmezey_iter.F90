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

subroutine PipekMezey_Iter(Functional,CMO,Ovlp,Thrs,ThrRot,ThrGrad,PA,nBas_per_Atom,nBas_Start,Name,nBasis,nOrb2Loc,nAtoms, &
                           nMxIter,Maximisation,Converged,Debug,Silent)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8 CMO(nBasis,*), Ovlp(nBasis,*)
real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
integer nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
character*(LENIN8) Name(nBasis)
logical Maximisation, Converged, Debug, Silent
real*8, allocatable :: RMat(:,:), PACol(:,:)

! Print iteration table header.
! -----------------------------

if (.not. Silent) then
  write(6,'(//,1X,A,/,1X,A)') '                                                        CPU       Wall', &
                              'nIter       Functional P        Delta     Gradient     (sec)     (sec) %Screen'
end if

! Initialization (iteration 0).
! -----------------------------

if (.not. Silent) call CWTime(C1,W1)
nIter = 0
call mma_Allocate(RMat,nOrb2Loc,nOrb2Loc,Label='RMat')
call GenerateP(Ovlp,CMO,Name,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Debug)
call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug)
OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional
if (.not. Silent) then
  call CWTime(C2,W2)
  TimC = C2-C1
  TimW = W2-W1
  write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,Zero
end if

! Iterations.
! -----------

call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')
Converged = .false.
do while ((nIter < nMxIter) .and. (.not. Converged))
  if (.not. Silent) call CWTime(C1,W1)
  call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,Maximisation,nOrb2Loc,Name,nBas_per_Atom,nBas_Start,ThrRot,PctSkp,Debug)
  call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
  call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,Debug)
  nIter = nIter+1
  Delta = Functional-OldFunctional
  OldFunctional = Functional
  if (.not. Silent) then
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,PctSkp
  end if
  Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)
end do
call mma_Deallocate(PACol)
call mma_Deallocate(RMat)

! Print convergence message.
! --------------------------

if (.not. Silent) then
  if (.not. Converged) then
    write(6,'(/,A,I4,A)') 'No convergence after',nIter,' iterations.'
  else
    write(6,'(/,A,I4,A)') 'Convergence after',nIter,' iterations.'
    write(6,*)
    write(6,'(A,I8)') 'Number of localised orbitals  : ',nOrb2loc
    write(6,'(A,1P,D20.10)') 'Value of P before localisation: ',FirstFunctional
    write(6,'(A,1P,D20.10)') 'Value of P after localisation : ',Functional
  end if
end if

end subroutine PipekMezey_Iter
