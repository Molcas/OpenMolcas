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

subroutine PipekMezey_Iter(Functional,CMO,Ovlp,PA,nBas_per_Atom,nBas_Start,BName,nBasis,nOrb2Loc,nAtoms,Converged)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.

use Molcas, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Pi
use Definitions, only: wp, iwp, u6
use Localisation_globals, only: Debug, Thrs,ThrGrad, Silent, nMxIter, OptMeth

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nBas_per_Atom(nAtoms), nBas_Start(nAtoms), nBasis, nOrb2Loc
real(kind=wp), intent(out) :: Functional, PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
character(len=LenIn+8), intent(in) :: BName(nBasis)
real(kind=wp), intent(in) :: Ovlp(nBasis,*)


logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: nIter, i,k, iBas, cnt, lSCR
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2, factor, ithrsh, DD, Thr
real(kind=wp), allocatable :: RMat(:,:), PACol(:,:),kappa(:,:),kappa_cnt(:,:),xkappa_cnt(:,:), &
                                GradientList(:,:,:), Hdiag(:,:), FunctionalList(:),&
                                unitary_mat(:,:), rotated_CMO(:,:), Ovlp_sqrt(:,:),  Ovlp_aux(:,:), SCR(:)
logical(kind=iwp), parameter :: printmore = .false., debug_exp = .false., debug_lowdin = .false.
real(kind=wp), parameter :: thrsh_taylor = 1.0e-16_wp, alpha = 0.3
real(kind=wp), External :: DDot_

! Initialization (iteration 0).
! -----------------------------

if (.not. Silent) call CWTime(C1,W1)

call mma_Allocate(RMat,nOrb2Loc,nOrb2Loc,Label='RMat')
call mma_Allocate(GradientList,nOrb2Loc,nOrb2Loc,nMxIter,Label='GradientList')  !nMxIter=300, maybe we can make it smaller
call mma_Allocate(FunctionalList,nMxIter,Label='FunctionalList')
call mma_Allocate(Hdiag,nOrb2Loc,nOrb2Loc,Label='Hdiag')
call mma_Allocate(unitary_mat,nOrb2Loc,nOrb2Loc,Label='unitary_mat')

nIter = 0
FunctionalList(:)=0

!if (lowdin) then
    call mma_allocate(Ovlp_sqrt, nBasis, nBasis,Label = "S^{1/2}")
    call mma_allocate(Ovlp_aux, nBasis, nBasis,Label = "S^{-1/2}")
    lSCR = 2*nBasis**2+nBasis*(nBasis+1)/2
    call mma_allocate(SCR,lSCR, Label = "SCR")
    if (debug_lowdin) then; call RecPrt("S before taking the sqrt",' ',Ovlp,nBasis,nBasis); end if
    call SQRTMT(Ovlp,nBasis,1,Ovlp_sqrt,Ovlp_aux,SCR)
    if (debug_lowdin) then
        call RecPrt("S^{1/2}",' ',Ovlp_sqrt,nBasis, nBasis)
        call RecPrt("S after taking the sqrt",' ',Ovlp,nBasis, nBasis)
        Ovlp_aux(:,:) = Zero ! i want to reuse it
        call dgemm_('N','N',nBasis,nBasis,nBasis,One,Ovlp_sqrt,nBasis,Ovlp_sqrt,nBasis,Zero,&
                    Ovlp_aux,nBasis)
        call RecPrt("S^{1/2}*S^{1/2}",' ',Ovlp_aux,nBasis,nBasis) ! should be same as S
    end if
    call mma_deallocate(Ovlp_aux)
!end if
call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)

if (.not. Silent) then
    write(u6,"(/A)") "MO extension before localisation:"
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)
else
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
end if

FunctionalList(1)=Functional

call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat, GradientList(:,:,1), Hdiag(:,:))

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
call mma_Allocate(rotated_cmo,nBasis,nOrb2Loc,Label='rotated_cmo') !this contains only the orbitals that are to be localized, while CMO always contains all
call mma_Allocate(kappa,nOrb2Loc,nOrb2Loc,Label='kappa')
call mma_Allocate(kappa_cnt,nOrb2Loc,nOrb2Loc,Label='kappa_cnt') != kappa^cnt
call mma_Allocate(xkappa_cnt,nOrb2Loc,nOrb2Loc,Label='xkappa_cnt') !saves the previous kappa_cnt

Converged = .false.

if (printmore) then
    write(u6,'(/,A)') '               nIter:  Functional:'
    write(u6,*) nIter,FunctionalList(nIter+1)
end if

do while ((nIter < nMxIter) .and. (.not. Converged) .and. (Functionallist(niter+1)<10*norb2loc))
    if (.not. Silent) call CWTime(C1,W1)

    !choose between optimization methods
    if (OptMeth == 1) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2x2 rotations
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,PctSkp)
    else if (OptMeth == 2 .or. OptMeth == 3) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! NxN rotations
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! define the transformation matrix

        kappa(:,:) = Zero
        kappa_cnt(:,:) = Zero
        xkappa_cnt(:,:) = Zero

        if (OptMeth == 2) then
            kappa(:,:) = -GradientList(:,:,nIter+1)/Hdiag(:,:)
        else if (OptMeth == 3) then
            kappa(:,:) = alpha*GradientList(:,:,nIter+1)
        end if
        DD=Sqrt(DDot_(nOrb2Loc**2,Kappa,1,Kappa,1))
        Thr= 0.5E0_wp * Pi
        If (DD>=Thr)Then
!           Write(6,*) 'Rescale Kappa(:,:)'
            Kappa(:,:) = (Thr/DD)*Kappa(:,:)
        End If

        kappa_cnt(:,:) = kappa !kappa^cnt = kappa since cnt=1
        xkappa_cnt(:,:) = kappa_cnt

        unitary_mat(:,:) = Zero
        call unitmat(unitary_mat,nOrb2Loc)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! analogous to exp_series in scf
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        cnt = 1
        factor = One
        ithrsh = 2.0e-16_wp

        unitary_mat(:,:) =  unitary_mat(:,:) - kappa(:,:)

        if (debug_exp) then
            write(u6,*) 'Taylor expansion: n=1; current iteration = ', nIter
            call RecPrt('unitary_mat = I - kappa^1',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
            call RecPrt('kappa',' ',kappa(:,:), nOrb2Loc, nOrb2Loc)
            write(u6,*) 'Taylor expansion: more terms'
        end if


        do while (ithrsh > thrsh_taylor)

            !the number of the term = the exponent for kappa in that term
            cnt = cnt+1

            !the faculty value that the matrix will be divided by
            factor = factor*DBLE(cnt)

            !calculate the cnt'th exponent of the kappa matrix
            ! this works by multiplying the matrix from the previous term (kappa^cnt) by the
            ! initial kappa matrix (just kappa^1)
            ! C <= alpha*A*B + beta*C
            ! kappa_cnt <= 1*kappa_cnt*kappa + 0*kappa_cnt
            call dgemm_('N','N',nOrb2Loc,nOrb2Loc,nOrb2Loc,(One/DBLE(cnt)),xkappa_cnt,nOrb2Loc,kappa,nOrb2Loc,Zero,&
                        kappa_cnt,norb2Loc)
            xkappa_cnt(:,:) = kappa_cnt

            ! differentiation of odd and even cases, because this expands exp(-kappa)
            ! all terms starting at n=2
            if (mod(cnt,2) == 0) then
                unitary_mat(:,:) =  unitary_mat + kappa_cnt(:,:)
                if (debug_exp) then
                    write(u6,'(A,F10.1,A,I2,A,ES12.4)') 'term: + 1/',factor,' * kappa^',cnt, &
                    ', current ithrsh = ', ithrsh
                end if
            else
                unitary_mat(:,:) =  unitary_mat - kappa_cnt(:,:)
                if (debug_exp) then
                    write(u6,'(A,F10.1,A,I2,A,ES12.4)') 'term: - 1/',factor,' * kappa^',cnt, &
                    ', current ithrsh = ', ithrsh
                end if
            end if

            ithrsh = maxval(abs(Kappa_Cnt(:,:))/(abs(unitary_mat)+thrsh_taylor))

            if (debug_exp) then
                write(u6,'(A,F10.1,A,I2,A,ES12.4)') 'term: + 1/',factor,' * kappa^',cnt, &
                    ', current ithrsh = ', ithrsh
                call RecPrt('kappa^cnt',' ',kappa_cnt(:,:), nOrb2Loc, nOrb2Loc)
                call RecPrt('unitary_mat',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
            end if
        end do

        if (debug) then
            write(u6,"(//A)") "rotating the orbitals with:"
            call RecPrt('kappa',' ',kappa(:,:), nOrb2Loc, nOrb2Loc)
            call RecPrt('unitary transformation matrix (exp(-kappa))',' ',unitary_mat(:,:), nOrb2Loc, nOrb2Loc)
            write(u6,"(//A,I5,/A)") "Starting next iteration ",niter+1,"--------------"
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! transform the orbitals
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rotated_CMO(:,:) = Zero
        !rotated_CMO(:,:) = CMO(:,:nOrb2Loc)

        !call dgemm_('N','N', nBasis, nOrb2Loc,nOrb2Loc,One, CMO(:,:nOrb2Loc), nOrb2Loc, unitary_mat, nOrb2Loc,&
        !            Zero,rotated_CMO,nOrb2Loc)
        do iBas = 1, nBasis
            do k = 1,nOrb2Loc
                do i = 1,nOrb2Loc
                    rotated_CMO(iBas,k) = rotated_CMO(iBas,k) + CMO(iBas,i) * unitary_mat(i,k)
                end do
            end do
        end do

        !reset CMO to be updated
        CMO(:,:) = rotated_CMO(:,:)
        call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)
    end if

    nIter = nIter+1

    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
    FunctionalList(nIter+1)=Functional !first entry is from before first iteration
    if (printmore) then
!       write(u6,'(/,A)') '               nIter:  Functional:'
        write(u6,*) nIter,FunctionalList(nIter+1)
    end if

    !calculates nxn gradient matrix for the current iteration and adds it to the List
    call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,RMat,GradientList(:,:,nIter+1), Hdiag(:,:))

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

call mma_Deallocate(kappa)
call mma_Deallocate(rotated_CMO)
call mma_Deallocate(kappa_cnt)
call mma_Deallocate(xkappa_cnt)
call mma_Deallocate(unitary_mat)
call mma_Deallocate(PACol)
call mma_Deallocate(RMat)
call mma_Deallocate(GradientList)
call mma_Deallocate(FunctionalList)
call mma_Deallocate(Hdiag)
!if (lowdin) then
    call mma_deallocate(Ovlp_sqrt)
    call mma_deallocate(SCR)
!end if

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
