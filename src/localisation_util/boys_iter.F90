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
subroutine Boys_Iter(Functional,CMO,Thrs,ThrRot,ThrGrad,ipLbl_AO,ipLbl,nBas,nOrb2Loc,nComp,nMxIter,Maximisation,Converged,Debug, &
                     Silent)
! Author: T.B. Pedersen
!
! Purpose: Boys localisation of orbitals.

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: Functional
real(kind=wp), intent(inout) :: CMO(*)
real(kind=wp), intent(in) :: Thrs, ThrRot, ThrGrad
integer(kind=iwp), intent(in) :: nComp, ipLbl_AO(nComp), ipLbl(nComp), nBas, nOrb2Loc, nMxIter
logical(kind=iwp), intent(in) :: Maximisation, Debug, Silent
logical(kind=iwp), intent(out) :: Converged
#include "WrkSpc.fh"
integer(kind=iwp) :: ipCol, ipRmat, lCol, lRmat, nIter
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2

! Print iteration table header.
! -----------------------------

if (.not. Silent) then
  write(u6,'(//,1X,A,/,1X,A)') '                                                        CPU       Wall', &
                               'nIter       Functional P        Delta     Gradient     (sec)     (sec) %Screen'
end if

! Initialization (iter 0).
! ------------------------

if (.not. Silent) call CWTime(C1,W1)
nIter = 0
Converged = .false.
lRmat = nOrb2Loc**2
call GetMem('Rmat','Allo','Real',ipRmat,lRmat)
call GenerateB(CMO,nBas,nOrb2Loc,ipLbl_AO,ipLbl,nComp,Debug)
call ComputeFuncB2(nOrb2Loc,ipLbl,nComp,Functional,Debug)
call GetGrad_Boys(nOrb2Loc,ipLbl,nComp,Work(ipRmat),GradNorm,Debug)
OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional
if (.not. Silent) then
  call CWTime(C2,W2)
  TimC = C2-C1
  TimW = W2-W1
  write(u6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,Zero
end if

! Iterations.
! -----------

lCol = 2*nOrb2Loc
call GetMem('Col','Allo','Real',ipCol,lCol)
do while ((nIter < nMxIter) .and. (.not. Converged))
  if (.not. Silent) call CWTime(C1,W1)
  call RotateOrbB(CMO,Work(ipCol),ipLbl,nComp,nBas,nOrb2Loc,Maximisation,ThrRot,PctSkp,Debug)
  call ComputeFuncB2(nOrb2Loc,ipLbl,nComp,Functional,Debug)
  call GetGrad_Boys(nOrb2Loc,ipLbl,nComp,Work(ipRmat),GradNorm,Debug)
  nIter = nIter+1
  Delta = Functional-OldFunctional
  OldFunctional = Functional
  if (.not. Silent) then
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(u6,'(1X,I5,1X,F18.8,2(1X,D12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,PctSkp
  end if
  Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)
end do
call GetMem('Col','Free','Real',ipCol,lCol)
call GetMem('Rmat','Free','Real',ipRmat,lRmat)

! Print convergence message.
! --------------------------

if (.not. Silent) then
  if (.not. Converged) then
    write(u6,'(/,A,I4,A)') 'No convergence after',nIter,' iterations.'
  else
    write(u6,'(/,A,I4,A)') 'Convergence after',nIter,' iterations.'
    write(u6,*)
    write(u6,'(A,1X,I4)') 'Number of localised orbitals  :',nOrb2Loc
    write(u6,'(A,1X,1P,D20.10)') 'Value of P before localisation:',FirstFunctional
    write(u6,'(A,1X,1P,D20.10)') 'Value of P after localisation :',Functional
  end if
end if

end subroutine Boys_Iter
