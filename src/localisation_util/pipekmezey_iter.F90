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

#define _DEBUGPRINT_
!#define _DEBUG2_
!#define _DEBUGLOWD_


subroutine PipekMezey_Iter(Functional,CMO,Ovlp,PA,nBas_per_Atom,nBas_Start,BName,nBasis,nOrb2Loc,nAtoms,Converged)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Pi
use Definitions, only: wp, iwp, u6
use Molcas, only: LenIn
use Localisation_globals, only: Thrs,ThrGrad, Silent, nMxIter, OptMeth, ChargeType, Loosen

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nBas_per_Atom(nAtoms), nBas_Start(nAtoms), nBasis, nOrb2Loc
real(kind=wp), intent(out) :: Functional, PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(in) :: Ovlp(nBasis,*)
character(len=LenIn+8), intent(in) :: BName(nBasis)
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: nIter, lSCR, fsdim
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2, DD, Thr,ang
real(kind=wp), allocatable :: PACol(:,:), GradientList(:,:), Functionallist(:), Hdiag(:,:), Ovlp_aux(:,:), &
                              SCR(:), Ovlp_sqrt(:,:),displacements(:,:),Gradient(:,:),dq(:),&
                              kappa(:,:),kappa_cnt(:,:),xkappa_cnt(:,:), unitary_mat(:,:), rotated_CMO(:,:),hdiagvec(:),&
                              Prev(:),Disp(:)
real(kind=wp), parameter :: alpha = 0.3
real(kind=wp), External :: DDot_
real(kind=wp) :: CtS(nOrb2Loc,nBasis),CtSC(nOrb2Loc,nOrb2Loc)

!S-GEK
real(kind=wp) :: dqdq
logical(kind=iwp) :: SORange
character(len=6):: UpMeth

! Initialization (iteration 0).
! -----------------------------

if (.not. Silent) call CWTime(C1,W1)

#ifdef _DEBUG2_
write(u6,'(/A)') 'Check the orthonormality of the orbitals'
write(u6,*) '========================================'
call dgemm_('T','N',nOrb2Loc, nBasis, nBasis,One, CMO, nBasis,Ovlp, nBasis,Zero, CtS, nOrb2Loc)
call dgemm_('N','N',nOrb2Loc, nOrb2Loc, nBasis,One,CtS, nOrb2Loc,CMO, nBasis,Zero,CtSC, nOrb2Loc)
call RecPrt("C^T*S*C =",' ',CtSC,nOrb2Loc, nOrb2Loc)
#endif

! to allow property printing later
call Put_cArray('Relax Method','LOCALIS ',8)

! allocating matrices for NxN optimizations
! ---------------------------------------------------------------------------------------------------
if (OptMeth == 2 .or. OptMeth == 3 .or. OptMeth == 4 .or. OptMeth == 5) then

    fsdim = nOrb2Loc*(nOrb2Loc-1)/2

    call mma_Allocate(kappa,nOrb2Loc,nOrb2Loc,Label='kappa')
    call mma_Allocate(Gradient,nOrb2Loc,nOrb2Loc,Label='Gradient')
    call mma_Allocate(Hdiag,nOrb2Loc,nOrb2Loc,Label='Hdiag')

    call mma_Allocate(Hdiagvec,fsdim,Label='Hdiagvec')
    call mma_Allocate(displacements,fsdim,nMxIter,Label='displacements')  ! kappa matrices
    call mma_Allocate(dq,fsdim,Label='dq')  ! GEK suggestion for kappa
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

#   ifdef _DEBUGLOWD_
        call RecPrt("S before taking the sqrt",' ',Ovlp,nBasis,nBasis)
#   endif

    call mma_allocate(SCR,lSCR, Label = "SCR")
    call SQRTMT(Ovlp,nBasis,1,Ovlp_sqrt,Ovlp_aux,SCR)
    call mma_Deallocate(SCR)

#   ifdef _DEBUGLOWD_
        call RecPrt("S^{1/2}",' ',Ovlp_sqrt,nBasis, nBasis)
        call RecPrt("S after taking the sqrt",' ',Ovlp,nBasis, nBasis)
        Ovlp_aux(:,:) = Zero ! i want to reuse it
        call dgemm_('N','N',nBasis,nBasis,nBasis,One,Ovlp_sqrt,nBasis,Ovlp_sqrt,nBasis,Zero,&
                    Ovlp_aux,nBasis)
        call RecPrt("S^{1/2}*S^{1/2}",' ',Ovlp_aux,nBasis,nBasis) ! should be same as S
#   endif

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
if (OptMeth == 2 .or. OptMeth == 3 .or. OptMeth == 4 .or. OptMeth == 5) then
    call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm, Gradient(:,:), Hdiagvec(:))
    call upper_triag2vec(Gradient(:,:),nOrb2Loc,GradientList(:,1),fsdim)
#   ifdef _DEBUGPRINT_
    call RecPrt("initial gradient"," ",Gradient,nOrb2Loc,nOrb2Loc)
    call RecPrt("initial hessian"," ",Hdiagvec(:),fsdim,1)
#   endif
    FunctionalList(1) = Functional
end if

! Print iteration table header.
! ---------------------------------------------------------------------------------------------------
OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional

UpMeth=" -  - "
if (.not. Silent) then
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(u6,'(//,1X,A,/,1X,A)') '                                                                   CPU       Wall', &
                                 'nIter       Functional P        Delta     Gradient   Microiter   (sec)     (sec) %Screen'
    write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),3X,A6,1X,2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,UpMeth,&
                                                    TimC,TimW,Zero
end if


! Iterations.
! ---------------------------------------------------------------------------------------------------
call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')
Converged = .false.

do while ((nIter < nMxIter) .and. (.not. Converged))
    if (.not. Silent) call CWTime(C1,W1)

    nIter = nIter+1

    !choose between optimization methods
    select case (OptMeth)

    case (1) ! Jacobi Sweeps (2x2 rotations)

        call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,PctSkp)
        call GetGradnorm_PM(nAtoms,nOrb2Loc,PA,GradNorm)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)

    case (2,3,4,5) ! Employing NxN rotations

        kappa(:,:) = Zero

        select case(OptMeth) !different NxN rot based methods

        case (2) ! Newton Raphson
            call vec2upper_triag(Hdiag(:,:),nOrb2Loc,Hdiagvec(:),fsdim,.false.)
            kappa(:,:) = -Gradient(:,:)/Hdiag(:,:)

        case (3) ! Gradient Ascent (no line search yet)
            kappa(:,:) = alpha*Gradient(:,:)

        case (4,5) ! (S)-GEK

            ! compute standard newton raphson step
            call vec2upper_triag(Hdiag(:,:),nOrb2Loc,Hdiagvec(:),fsdim,.false.)
            kappa(:,:) = -Gradient(:,:)/Hdiag(:,:)
            call upper_triag2vec(kappa(:,:),nOrb2Loc,displacements(:,nIter+1),fsdim)

#           ifdef _DEBUGPRINT_
                write(u6,"(A,I3,A)") "kappa_",nIter,"="
                call RecPrt('-g/hdiag (NR step) as matrix',' ',kappa(:,:),nOrb2Loc,nOrb2Loc)
                call RecPrt('-g/hdiag (NR step) as vector',' ',displacements(:,nIter+1),fsdim,1)
#           endif

            !if (nIter == 250) then
            !    call get_intermediate_molden(CMO,nBasis,nOrb2Loc)
            !end if

            ! start GEK only from iteration x
            if (nIter > 1) then

                SORange = .true.

                select case(OptMeth)

                case (4) ! Full space GEK

                    call S_GEK_localisation(nIter,Functionallist(:),-GradientList(:,:),displacements(:,:),-hdiagvec(:),fsdim,&
                                            dqdq,displacements(:,nIter+1),UpMeth,'fullspace',SORange)

                case (5) ! subspace GEK

                    write(u6,*) "building the subspace"
                    call S_GEK_localisation(nIter,Functionallist(:),-GradientList(:,:),displacements(:,:),-hdiagvec(:),fsdim,&
                                        dqdq,displacements(:,nIter+1),UpMeth,'subspace ',SORange)
                end select !(s)-GEK


                ! undershoot mitigation
                if (Loosen%Step > One) then
                    call mma_allocate(Prev,fsdim,Label='Prev')
                    call mma_allocate(Disp,fsdim,Label='Disp')

                    Prev(:) = displacements(:,nIter)
                    Disp(:) = displacements(:,nIter+1)

                    dqdq = DDot_(fsdim,Disp,1,Disp,1)*DDot_(fsdim,Prev,1,Prev,1)
                    ang = DDot_(fsdim,Prev,1,Disp,1)/sqrt(dqdq)
                    if (ang < Loosen%Thrs2) then
                        Loosen%Factor = One
                    else if (ang > Loosen%Thrs) then
                        Loosen%Factor = Loosen%Factor*Loosen%Step
                    end if

#                   ifdef _DEBUGPRINT_
                        call RecPrt('Disp',' ',Disp,fsdim,1)
                        call RecPrt('Prev',' ',Prev,fsdim,1)
                        write(u6,*) "angle(Disp,Prev) = cos^-1(",ang,")"
                        write(u6,*) "Loosen%Factor    =", Loosen%Factor
                        write(u6,*) "Loosen%Step    =", Loosen%Step
#                   endif

                    call mma_Deallocate(Prev)
                    call mma_Deallocate(Disp)

                end if

                ! transform GEK disp vec to matrix
                call vec2upper_triag(kappa(:,:),nOrb2Loc,displacements(:,nIter+1),fsdim,.true.)
#               ifdef _DEBUGPRINT_
                call RecPrt('(GEK step)',' ',displacements(:,nIter+1),fsdim,1)
#               endif


            end if
        end select ! different NxN rotations
        ! ---------------------------------------------------------------------------------------------------

        DD=Sqrt(DDot_(nOrb2Loc**2,Kappa,1,Kappa,1))
        !Thr= 0.5E0_wp * Pi
        Thr= Pi
        If (DD>=Thr)Then
#           ifdef _DEBUGPRINT_
            Write(6,*) 'Rescale Kappa(:,:)'
#           endif
            Kappa(:,:) = (Thr/DD)*Kappa(:,:)
        End If


#       ifdef _DEBUGPRINT_
            call RecPrt('displacement taken (kappa mat)',' ',kappa(:,:),nOrb2Loc,nOrb2Loc)
#       endif

        call RotateNxN(CMO,kappa,nOrb2Loc,nBasis,kappa_cnt,xkappa_cnt,unitary_mat,rotated_CMO)

#       ifdef _DEBUGPRINT_
            write(u6,*) "=================================================================="
            write(u6,*) "               ORBITALS HAVE BEEN ROTATED"
            write(u6,*) "=================================================================="
#       endif

        call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Gradient(:,:), Hdiagvec(:)) ! gets the new gradient

#       ifdef _DEBUGPRINT_
            write(u6,*) "               NEW GRADIENT & NEW HESSIAN DIAGONAL:               "
            write(u6,*) "Functional:",Functional
            call RecPrt('Gradient',' ',Gradient(:,:),nOrb2Loc,nOrb2Loc)
            call RecPrt('Hdiag',' ',Hdiagvec(:),fsdim,1)
#       endif

        call upper_triag2vec(Gradient(:,:),nOrb2Loc,GradientList(:,nIter+1),fsdim)

        call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
        FunctionalList(nIter+1)=Functional !first entry is from before first iteration

    end select ! 2x2 or NxN rotations

    !check if converged
    ! ---------------------------------------------------------------------------------------------------
    Delta = Functional-OldFunctional
    OldFunctional = Functional
    if (.not. Silent) then
        call CWTime(C2,W2)
        TimC = C2-C1
        TimW = W2-W1
        write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),3X,A6,1X,2(1X,F9.1),1X,F7.2)') nIter,Functional,Delta,GradNorm,UpMeth,&
                                                                                TimC,TimW,PctSkp
    end if
    Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)

    !this is just to see the orbitals (REMOVE LATER)
    if (nIter == nMxIter-2) then
        Converged = .true.
    end if
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

!call Prpt()
#ifdef _DEBUGPRINT_
write(u6,'(/A)') 'Check the orthonormality of the orbitals'
write(u6,*) '========================================'
call dgemm_('T','N',nOrb2Loc, nBasis, nBasis,One, CMO, nBasis,Ovlp, nBasis,Zero, CtS, nOrb2Loc)
call dgemm_('N','N',nOrb2Loc, nOrb2Loc, nBasis,One,CtS, nOrb2Loc,CMO, nBasis,Zero,CtSC, nOrb2Loc)
call RecPrt("C^T*S*C =",' ',CtSC,nOrb2Loc, nOrb2Loc)
#endif


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
    call mma_Deallocate(Hdiagvec)
    call mma_Deallocate(dq)
end if

! deallocate other matrices
call mma_Deallocate(PACol)
call mma_Deallocate(Ovlp_sqrt)

end subroutine PipekMezey_Iter
