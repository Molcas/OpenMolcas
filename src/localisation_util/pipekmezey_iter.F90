!**********************************************************************
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

subroutine PipekMezey_Iter(Functional,CMO,Ovlp,PA,nBas_per_Atom,nBas_Start,BName,nBasis,nOrb2Loc,nAtoms,Converged)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Pi
use Definitions, only: wp, iwp, u6
use Molcas, only: LenIn8
use Localisation_globals, only: Debug, Thrs,ThrGrad, Silent, nMxIter, OptMeth, ChargeType

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nBas_per_Atom(nAtoms), nBas_Start(nAtoms), nBasis, nOrb2Loc
real(kind=wp), intent(out) :: Functional, PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(in) :: Ovlp(nBasis,*)
character(len=LenIn8), intent(in) :: BName(nBasis)
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: nIter, i,k, iBas, lSCR
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2
real(kind=wp), allocatable :: PACol(:,:), GradientList(:,:,:), Functionallist(:), Hdiag(:,:), Ovlp_aux(:,:), &
                              SCR(:), Ovlp_sqrt(:,:)
logical(kind=iwp), parameter :: printmore = .false., debug_lowdin = .false.

! Initialization (iteration 0).
! -----------------------------

if (.not. Silent) call CWTime(C1,W1)
call mma_Allocate(GradientList,nOrb2Loc,nOrb2Loc,nMxIter,Label='GradientList')  !nMxIter=300, maybe we can make it smaller
call mma_Allocate(FunctionalList,nMxIter,Label='FunctionalList')
call mma_Allocate(Hdiag,nOrb2Loc,nOrb2Loc,Label='Hdiag')
nIter = 0
FunctionalList(:)=0

call mma_allocate(Ovlp_sqrt, nBasis, nBasis,Label = "S^{1/2}")

if (ChargeType ==2) then !Lowdin
    call mma_allocate(Ovlp_aux, nBasis, nBasis,Label = "S^{-1/2}")
    lSCR = 2*nBasis**2+nBasis*(nBasis+1)/2
    if (debug_lowdin) then; call RecPrt("S before taking the sqrt",' ',Ovlp,nBasis,nBasis); end if
    call mma_allocate(SCR,lSCR, Label = "SCR")
    call SQRTMT(Ovlp,nBasis,1,Ovlp_sqrt,Ovlp_aux,SCR)
    call mma_Deallocate(SCR)
    if (debug_lowdin) then
        call RecPrt("S^{1/2}",' ',Ovlp_sqrt,nBasis, nBasis)
        call RecPrt("S after taking the sqrt",' ',Ovlp,nBasis, nBasis)
        Ovlp_aux(:,:) = Zero ! i want to reuse it
        call dgemm_('N','N',nBasis,nBasis,nBasis,One,Ovlp_sqrt,nBasis,Ovlp_sqrt,nBasis,Zero,&
                    Ovlp_aux,nBasis)
        call RecPrt("S^{1/2}*S^{1/2}",' ',Ovlp_aux,nBasis,nBasis) ! should be same as S
    end if
    call mma_deallocate(Ovlp_aux)
end if

call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)

if (.not. Silent) then
    write(u6,"(/A)") "MO extension before localisation:"
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)
else
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
end if

FunctionalList(1)=Functional

call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm, GradientList(:,:,1), Hdiag(:,:))

OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional

! Print iteration table header.
! -----------------------------

if (.not. Silent) then
  write(u6,'(//,1X,A,/,1X,A)') '                                                        CPU       Wall', &
                               'nIter       Functional P        Delta     Gradient     (sec)     (sec) %Screen'
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,Zero
end if

! Iterations.
! -----------

call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')
Converged = .false.

if (printmore) then
    write(u6,'(/,A)') '               nIter:  Functional:'
    write(u6,*) nIter,FunctionalList(nIter+1)
end if

do while ((nIter < nMxIter) .and. (.not. Converged) .and. (Functionallist(niter+1)<10*norb2loc))
    if (.not. Silent) call CWTime(C1,W1)
    !choose between optimization methods
    if (OptMeth == 1) then
        ! 2x2 rotations: Jacobi Sweeps
        call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,PctSkp)
    else if (OptMeth == 2 .or. OptMeth == 3) then
        ! NXN rotations: Gradient Ascent or Newton Raphson
        call RotateNxN(CMO,Ovlp,nOrb2Loc,nBasis,Ovlp_sqrt(:,:),GradientList(:,:,nIter+1),Hdiag(:,:),BName,nAtoms,&
                       nBas_per_Atom,nBas_Start,PA(:,:,:))
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,GradientList(:,:,nIter+2), Hdiag(:,:))
   end if

    nIter = nIter+1

    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
    FunctionalList(nIter+1)=Functional !first entry is from before first iteration
    if (printmore) then
!       write(u6,'(/,A)') '               nIter:  Functional:'
        write(u6,*) nIter,FunctionalList(nIter+1)
    end if

    !calculates nxn gradient matrix for the current iteration and adds it to the List
    call GetGradnorm_PM(nAtoms,nOrb2Loc,PA,GradNorm)

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

! print info about each localized MO
if (.not. Silent) then
    write(u6,"(/A)") "MO extension after localisation:"
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)
end if

call mma_Deallocate(PACol)
call mma_Deallocate(GradientList)
call mma_Deallocate(FunctionalList)
call mma_Deallocate(Hdiag)
call mma_Deallocate(Ovlp_sqrt)

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

