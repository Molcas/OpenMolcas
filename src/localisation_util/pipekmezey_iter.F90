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
!               2026, Lila Zapp (opt methods & loewdin framework)      *
!***********************************************************************

subroutine PipekMezey_Iter(Functional,CMO,Ovlp,PA,nBas_per_Atom,nBas_Start,BName,nBasis,nOrb2Loc,nAtoms,Converged)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Pi
use Definitions, only: wp, iwp, u6
use Molcas, only: LenIn
use Localisation_globals, only: Thrs,ThrGrad, Silent, nMxIter, OptMeth, ChargeType

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nBas_per_Atom(nAtoms), nBas_Start(nAtoms), nBasis, nOrb2Loc
real(kind=wp), intent(out) :: Functional, PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(in) :: Ovlp(nBasis,*)
character(len=LenIn+8), intent(in) :: BName(nBasis)
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: nIter, lSCR
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2, DD, Thr
real(kind=wp), allocatable :: PACol(:,:), GradientList(:,:), Functionallist(:), Hdiag(:,:), Ovlp_aux(:,:), &
                              SCR(:), Ovlp_sqrt(:,:),displacements(:,:),Gradient(:,:),&
                              kappa(:,:),kappa_cnt(:,:),xkappa_cnt(:,:), unitary_mat(:,:), rotated_CMO(:,:)
logical(kind=iwp), parameter :: debug_lowdin = .false.
real(kind=wp), parameter :: alpha = 0.3
real(kind=wp), External :: DDot_

!for S-GEK
integer(kind=iwp) :: nDiis,iFirst,fsdim,i,j,k,l,nExplicit,mDiis
real(kind=wp) :: gg
real(kind=wp), allocatable :: q(:,:),g(:,:),Aux_a(:),Aux_b(:),e_diis(:,:),dq(:)
integer(kind=iwp), parameter :: nWindow = 5

! Initialization (iteration 0).
! -----------------------------

if (.not. Silent) call CWTime(C1,W1)


! allocating matrices for NxN optimizations
! ---------------------------------------------------------------------------------------------------
if (OptMeth == 2 .or. OptMeth == 3 .or. OptMeth == 4) then

    fsdim = nOrb2Loc*(nOrb2Loc-1)/2

    call mma_Allocate(kappa,nOrb2Loc,nOrb2Loc,Label='kappa')
    call mma_Allocate(Gradient,nOrb2Loc,nOrb2Loc,Label='Gradient')
    call mma_Allocate(Hdiag,nOrb2Loc,nOrb2Loc,Label='Hdiag')

    call mma_Allocate(displacements,fsdim,nMxIter,Label='displacements')  ! kappa matrices
    call mma_Allocate(GradientList,fsdim,nMxIter,Label='GradientList')
    call mma_Allocate(FunctionalList,nMxIter,Label='FunctionalList')
    displacements(:,:)=Zero
    GradientList(:,:)=Zero
    FunctionalList(:)=Zero

    call mma_Allocate(kappa_cnt,nOrb2Loc,nOrb2Loc,Label='kappa_cnt') != kappa^cnt
    call mma_Allocate(xkappa_cnt,nOrb2Loc,nOrb2Loc,Label='xkappa_cnt') !saves the previous kappa_cnt
    call mma_Allocate(unitary_mat,nOrb2Loc,nOrb2Loc,Label='unitary_mat')
    call mma_Allocate(rotated_cmo,nBasis,nOrb2Loc,Label='rotated_cmo')
end if
! ---------------------------------------------------------------------------------------------------

nIter = 0

call mma_allocate(Ovlp_sqrt, nBasis, nBasis,Label = "S^{1/2}")


! if the Loewdin charge framework is requested instead of Mulliken
! ---------------------------------------------------------------------------------------------------
if (ChargeType ==2) then
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
! ---------------------------------------------------------------------------------------------------

call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)


if (.not. Silent) then
    write(u6,"(/A)") "MO extension before localisation:"
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)
else
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
end if


! get initial gradient, hessian diagonal, add initial functional value to list
! ---------------------------------------------------------------------------------------------------
if (OptMeth == 2 .or. OptMeth == 3 .or. OptMeth == 4) then
    call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm, Gradient(:,:), Hdiag(:,:))
    call upper_triag2vec(Gradient(:,:),nOrb2Loc,GradientList(:,1),fsdim)
    call RecPrt("initial gradient"," ",Gradient,nOrb2Loc,nOrb2Loc)
    FunctionalList(1) = Functional
end if


! Print iteration table header.
! ---------------------------------------------------------------------------------------------------
OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional

if (.not. Silent) then
  write(u6,'(//,1X,A,/,1X,A)') '                                                        CPU       Wall', &
                               'nIter       Functional P        Delta     Gradient     (sec)     (sec) %Screen'
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,Zero
end if


! Iterations.
! ---------------------------------------------------------------------------------------------------
call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')
Converged = .false.

do while ((nIter < nMxIter) .and. (.not. Converged))
    if (.not. Silent) call CWTime(C1,W1)

    nIter = nIter+1

    !choose between optimization methods

    ! Jacobi Sweeps (2x2 rotations)
    ! ---------------------------------------------------------------------------------------------------
    if (OptMeth == 1) then
        call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,PctSkp)
        call GetGradnorm_PM(nAtoms,nOrb2Loc,PA,GradNorm)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)

    ! Employing NxN rotations
    ! ---------------------------------------------------------------------------------------------------
    else if (OptMeth == 2 .or. OptMeth == 3 .or. OptMeth == 4) then

        kappa(:,:) = Zero

        ! Newton Raphson
        ! ---------------------------------------------------------------------------------------------------
        if (OptMeth == 2) then
            kappa(:,:) = -Gradient(:,:)/Hdiag(:,:)


        ! Gradient Ascent (no line search yet)
        ! ---------------------------------------------------------------------------------------------------
        else if (OptMeth == 3) then
            kappa(:,:) = alpha*Gradient(:,:)


        ! S-GEK
        ! ---------------------------------------------------------------------------------------------------
        else if (OptMeth == 4) then ! S-GEK

            ! Pick up coordinates and gradients in full space
            ! -----------------------------------------------

            ! number of iterations used to build the subspace
            nDIIS = min(nIter,nWindow) !1 for first iteration; 2

            if (nDIIS == 1) then
                write(u6,*) 'Exit S-GEK Optimizer'
                kappa(:,:) = -Gradient(:,:)/Hdiag(:,:)
            else


            ! index of the first iteration to consider for the subspace
            iFirst = nIter-nDIIS+1 !1 for first iteration; 1

            call mma_Allocate(q,fsdim, nDiis,Label="q")
            call mma_Allocate(g,fsdim, nDiis,Label="g")

            call mma_Allocate(dq,fsdim,Label='dq')

            j = 0
            do i=iFirst,nIter
                j = i-iFirst+1
                !write(u6,*) 'i,j,iter=',i,j,nIter

                ! Coordinates
                q(:,j) = displacements(:,i)

                ! Gradients
                g(:,j) = GradientList(:,i)

            end do

            !change this later
            dq(:) = displacements(:,nIter)

            write(u6,*) 'nWindow =',nWindow
            write(u6,*) '  nDIIS =',nDIIS
            write(u6,*) '  nIter =',nIter
            call RecPrt("g(:,:)",' ',g,fsdim, nDiis)
            call RecPrt("q(:,:)",' ',q,fsdim, nDiis)
            call RecPrt("g(:,nDiis)",' ',g(:,nDiis),fsdim, 1)
            call RecPrt("dq(:)",' ',dq,fsdim, 1)

            ! select subspace basis vectors; construct normalized e_diis
            ! -----------------------------------------------------------

            !number of subspace basis vectors, potentially linear dependent => difference vecs of ndiis displacements and gradients +2 additional vecs (see below)
            nExplicit = 2*(nDIIS-1)+2

            call mma_allocate(e_diis,fsdim,nExplicit,Label='e_diis')

            call mma_allocate(Aux_a,fsdim,Label='Aux_a')
            call mma_allocate(Aux_b,fsdim,Label='Aux_b')

            j = 0
            do k=1,nDIIS-1
                !n-th column of e_diis
                j = j+1
                ! gradient difference vector
                Aux_a(:) = g(:,k+1)-g(:,k)
                !normalize
                e_diis(:,j) = Aux_a(:)/sqrt(DDot_(fsdim,Aux_a(:),1,Aux_a(:),1))

                !(n+1)-th column of e_diis
                j = j+1
                ! displacement difference vector
                Aux_a(:) = q(:,k+1)-q(:,k)
                Aux_b(:) = Aux_a(:)
                !normalize
                e_diis(:,j) = Aux_b(:)/sqrt(DDot_(fsdim,Aux_b(:),1,Aux_b(:),1))

            end do
            call mma_deallocate(Aux_b)

            ! Add some unit vectors corresponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
            j = j+1
            !current gradient
            Aux_a(:) = g(:,nDIIS)
            !normalize
            e_diis(:,j) = Aux_a(:)/sqrt(DDot_(fsdim,Aux_a(:),1,Aux_a(:),1))

            j = j+1
            !second order method's displacement suggestion
            Aux_a(:) = dq(:)
            !normalize
            e_diis(:,j) = Aux_a(:)/sqrt(DDot_(fsdim,Aux_a(:),1,Aux_a(:),1))
            call mma_deallocate(Aux_a)

            if (allocated(e_diis)) call RecPrt('e_diis(unorth)',' ',e_diis,fsdim,nExplicit)

            call mma_Deallocate(q)
            call mma_Deallocate(g)
            call mma_Deallocate(dq)


            ! orthogonalize e_diis; remove redundancies from linear dependences
            ! -----------------------------------------------------------------
            do l=1,2
                j = 1
                do i=2,nExplicit
                    do k=1,j
                        gg = DDot_(fsdim,e_diis(:,i),1,e_diis(:,k),1)
                        write(u6,*) 'i,k,gg=',i,k,gg
                        e_diis(:,i) = e_diis(:,i)-gg*e_diis(:,k)
                    end do
                    gg = DDot_(fsdim,e_diis(:,i),1,e_diis(:,i),1) ! renormalize
                    write(u6,*) 'j,i,gg=',j,i,gg

                    if (gg > 1.0e-17_wp) then   ! Skip vector if linear dependent.
                        j = j+1
                        e_diis(:,j) = e_diis(:,i)/sqrt(gg)
                    end if
                end do
            end do

            ! normally mDIIS=2*nDIIS, but it can happen that not all unit vectors are linear independent (mDIIS<=2*nDIIS).
            ! mDIIS is then the number of linear independent e_diis column vectors that span the subspace
            mDIIS = j


            write(u6,*) '    fsdim:',fsdim
            write(u6,*) 'nExplicit:',nExplicit
            write(u6,*) '    nIter:',nIter
            write(u6,*) '    nDIIS:',nDIIS
            write(u6,*) '    mDIIS:',mDIIS

            write(u6,*) 'Check the orthonormality'
            do i=1,mDIIS
                do j=1,i
                    write(u6,*) i,j,DDot_(fsdim,e_diis(:,i),1,e_diis(:,j),1)
                end do
                write(u6,*)
            end do
            if (allocated(e_diis)) call RecPrt('e_diis',' ',e_diis,fsdim,mDIIS)












            call mma_Deallocate(e_diis)


            !call S_GEK_localisation(nOrb2Loc,nDiis,kappa,GradientList(:,:),displacements(:,:))

            kappa(:,:) = -Gradient(:,:)/Hdiag(:,:)
            call upper_triag2vec(kappa(:,:),nOrb2Loc,displacements(:,nIter+1),fsdim)
            end if ! S-GEK with nIter > 1
        end if ! different NxN rotations
        ! ---------------------------------------------------------------------------------------------------

        DD=Sqrt(DDot_(nOrb2Loc**2,Kappa,1,Kappa,1))
        Thr= 0.5E0_wp * Pi
        If (DD>=Thr)Then
        !           Write(6,*) 'Rescale Kappa(:,:)'
            Kappa(:,:) = (Thr/DD)*Kappa(:,:)
        End If

        call RotateNxN(CMO,kappa,nOrb2Loc,nBasis,kappa_cnt,xkappa_cnt,unitary_mat,rotated_CMO)
        call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Gradient(:,:), Hdiag(:,:)) ! gets the new gradient
        call upper_triag2vec(Gradient(:,:),nOrb2Loc,GradientList(:,nIter+1),fsdim)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
        FunctionalList(nIter+1)=Functional !first entry is from before first iteration

    end if ! different opt methods

    !check if converged
    ! ---------------------------------------------------------------------------------------------------
    Delta = Functional-OldFunctional
    OldFunctional = Functional
    if (.not. Silent) then
        call CWTime(C2,W2)
        TimC = C2-C1
        TimW = W2-W1
        write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,TimC,TimW,PctSkp
    end if
    Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)
end do !Iterations


! print info about each localized MO
! ---------------------------------------------------------------------------------------------------
if (.not. Silent) then
    write(u6,"(/A)") "MO extension after localisation:"
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)
end if


! Print convergence message.
! ---------------------------------------------------------------------------------------------------
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

! deallocate matrices used for NxN optimizations
! ---------------------------------------------------------------------------------------------------
if (OptMeth == 2 .or. OptMeth == 3 .or. OptMeth == 4) then
    call mma_Deallocate(Gradient)
    call mma_Deallocate(Hdiag)
    call mma_Deallocate(kappa)

    call mma_Deallocate(kappa_cnt)
    call mma_Deallocate(xkappa_cnt)
    call mma_Deallocate(unitary_mat)
    call mma_Deallocate(rotated_CMO)

    call mma_Deallocate(FunctionalList)
    call mma_Deallocate(GradientList)
    call mma_Deallocate(displacements)
end if

! deallocate other matrices
call mma_Deallocate(PACol)
call mma_Deallocate(Ovlp_sqrt)

end subroutine PipekMezey_Iter
