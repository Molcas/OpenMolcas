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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************
subroutine EdmistonRuedenberg_Iter(Functional,CMO,Thrs,ThrRot,ThrGrad,nBasis,nOrb2Loc,nMxIter,Maximisation,Converged,Debug,Silent)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: ER localisation of orbitals.
!
! The optimization algorithm is the "generalized eta step" of
! Subotnik, Shao, Liang, and Head-Gordon, JCP 121, 9220 (2004).
!
! Redundant arguments that might be used at a later stage:
!   ThrRot [might be used for DIIS]
!   Maximisation [might be used for Jacobi sweeps]
!
! Note that two-electron integrals (Cholesky decomposed) must be
! available and appropriately set up when calling this routine.

implicit real*8(a-h,o-z)
real*8 CMO(nBasis,nOrb2Loc)
logical Maximisation, Converged, Debug, Silent
#include "WrkSpc.fh"

character*23 SecNam
parameter(SecNam='EdmistonRuedenberg_Iter')

logical Timing

if (Debug) then
  write(6,*) SecNam,'[debug]: Maximisation: ',Maximisation
  write(6,*) SecNam,'[debug]: ThrRot      : ',ThrRot
end if

! Print iteration table header.
! -----------------------------

if (.not. Silent) then
  write(6,'(//,1X,A,/,1X,A)') '                                                        CPU       Wall', &
                              'nIter      Functional ER        Delta     Gradient     (sec)     (sec)'
end if

! Initialization.
! ---------------

Converged = .false.
Timing = Debug

lRmat = nOrb2Loc**2
call GetMem('Rmat','Allo','Real',ipRmat,lRmat)

! Iteration 0.
! ------------

if (.not. Silent) call CWTime(C1,W1)
nIter = 0
Functional = 0.0d0
call GetGrad_ER(Functional,GradNorm,Work(ipRmat),CMO,nBasis,nOrb2Loc,Timing)
OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional
if (.not. Silent) then
  call CWTime(C2,W2)
  TimC = C2-C1
  TimW = W2-W1
  write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1))') nIter,Functional,Delta,GradNorm,TimC,TimW
end if

! Iterations.
! -----------

do while ((nIter < nMxIter) .and. (.not. Converged))
  if (.not. Silent) call CWTime(C1,W1)
  call RotateOrb_ER(Work(ipRmat),CMO,nBasis,nOrb2Loc,Debug)
  call GetGrad_ER(Functional,GradNorm,Work(ipRmat),CMO,nBasis,nOrb2Loc,Timing)
  nIter = nIter+1
  Delta = Functional-OldFunctional
  OldFunctional = Functional
  if (.not. Silent) then
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1))') nIter,Functional,Delta,GradNorm,TimC,TimW
  end if
  Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)
end do

! Print convergence message.
! --------------------------

if (.not. Silent) then
  if (.not. Converged) then
    write(6,'(/,A,I4,A)') 'No convergence after',nIter,' iterations.'
  else
    write(6,'(/,A,I4,A)') 'Convergence after',nIter,' iterations.'
    write(6,*)
    write(6,'(A,I8)') 'Number of localised orbitals  : ',nOrb2loc
    write(6,'(A,F12.8)') 'Value of P before localisation: ',FirstFunctional
    write(6,'(A,F12.8)') 'Value of P after localisation : ',Functional
  end if
end if

! Finalization.
! -------------

call GetMem('Rmat','Free','Real',ipRmat,lRmat)

end subroutine EdmistonRuedenberg_Iter
