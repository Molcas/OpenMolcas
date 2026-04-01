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

!#define _DEBUGLISTS_
!#define _DEBUG2_
!#define _DEBUGPRINT_
!#define _DEBUGLOWD_
!#define _GETMOLDEN_

subroutine PipekMezey_Iter(Functional,CMO,Ovlp,PA,nBas_per_Atom,nBas_Start,BName,nBasis,nOrb2Loc,nAtoms,Converged)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Pi
use Definitions, only: wp, iwp, u6
use Molcas, only: LenIn
use Localisation_globals, only: Thrs,ThrGrad, Silent, nMxIter, OptMeth, ChargeType, Loosen, FuncList, GradList, DispList,&
                                UmatList
#ifdef _GETMOLDEN_
use filesystem, only: getcwd_, mkdir_
#endif

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nBas_per_Atom(nAtoms), nBas_Start(nAtoms), nBasis, nOrb2Loc
real(kind=wp), intent(out) :: Functional, PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(in) :: Ovlp(nBasis,*)
character(len=LenIn+8), intent(in) :: BName(nBasis)
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: nIter, lSCR, fsdim,nDIIS
real(kind=wp) :: C1, C2, Delta, FirstFunctional, GradNorm, OldFunctional, PctSkp, TimC, TimW, W1, W2, ang
real(kind=wp), allocatable :: PACol(:,:), Hdiag(:,:), Ovlp_aux(:,:), &
                              SCR(:), Ovlp_sqrt(:,:),Gradient(:),&
                              kappa(:,:),kappa_cnt(:,:),xkappa_cnt(:,:), unitary_mat(:,:), rotated_CMO(:,:),hdiagvec(:),&
                              Prev(:),Disp(:),CMO_Ref(:,:)
real(kind=wp), parameter :: alpha = 0.3
real(kind=wp), External :: DDot_

! for S-GEK
integer(kind=iwp) :: maxel
real(kind=wp) :: dqdq,largest
logical(kind=iwp) :: SORange,GEKRange,ResetGEK
character(len=6):: UpMeth
logical(kind=iwp),parameter :: usmitigation = .false.
integer(kind=iwp) :: i,IterGEK,large_elements,mindp

real(kind=wp) :: DD,Thr
#ifdef _DEBUGPRINT_
real(kind=wp) :: CtS(nOrb2Loc,nBasis),CtSC(nOrb2Loc,nOrb2Loc)
#endif
#ifdef _GETMOLDEN_
character(len=1024) :: Sub, WorkDir, NewDir, SubmitDir, imfile
integer(kind=iwp) :: rc
character(len=8) :: fmt
character(len=4) :: x1
#endif


# ifdef _GETMOLDEN_

! preparations
fmt = '(I4.4)'

! locate scratch directory
call getcwd_(WorkDir)
write(u6,*) "WorkDir = ", trim(WorkDir)

! Create intermediate_molden directory that contains molden files of every iteration
Sub = "intermediate_molden"
call getenvf('MOLCAS_SUBMIT_DIR',SubmitDir)
NewDir = trim(SubmitDir)//'/'//Sub
call mkdir_(NewDir)

# endif


if (.not. Silent) call CWTime(C1,W1)

#ifdef _DEBUGPRINT_
write(u6,'(/A)') 'Check the orthonormality of the orbitals'
write(u6,*) '========================================'
call dgemm_('T','N',nOrb2Loc, nBasis, nBasis,One, CMO, nBasis,Ovlp, nBasis,Zero, CtS, nOrb2Loc)
call dgemm_('N','N',nOrb2Loc, nOrb2Loc, nBasis,One,CtS, nOrb2Loc,CMO, nBasis,Zero,CtSC, nOrb2Loc)
call RecPrt("C^T*S*C =",' ',CtSC,nOrb2Loc, nOrb2Loc)
#endif


! to allow property printing later
call Put_cArray('Relax Method','LOCALIS ',8)


! if the Loewdin charge framework is requested instead of Mulliken
! ---------------------------------------------------------------------------------------------------
call mma_allocate(Ovlp_sqrt, nBasis, nBasis,Label = "S^{1/2}")
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


! allocations

select case(OptMeth)

case (1)

    call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')

case(2,3,4,5)

    ! allocating matrices for NxN optimizations

    fsdim = nOrb2Loc*(nOrb2Loc-1)/2

    call mma_Allocate(kappa,nOrb2Loc,nOrb2Loc,Label='kappa')
    call mma_Allocate(Gradient,fsdim,Label='Gradient')
    call mma_Allocate(Hdiag,nOrb2Loc,nOrb2Loc,Label='Hdiag')

    call mma_Allocate(Hdiagvec,fsdim,Label='Hdiagvec')
    call mma_Allocate(DispList,fsdim,nMxIter,Label='DispList')  ! kappa matrices
    call mma_Allocate(UmatList,nOrb2Loc,nOrb2Loc,nMxIter,Label='UmatList')
    call mma_allocate(Disp,fsdim,Label='Disp')
    call mma_Allocate(GradList,fsdim,nMxIter,Label='GradList')
    call mma_Allocate(FuncList,nMxIter,Label='FuncList')
    DispList(:,:)=Zero
    UmatList(:,:,:)=Zero
    GradList(:,:)=Zero
    FuncList(:)=Zero
    Kappa(:,:)=Zero
    call mma_Allocate(kappa_cnt,nOrb2Loc,nOrb2Loc,Label='kappa_cnt') != kappa^cnt
    call mma_Allocate(xkappa_cnt,nOrb2Loc,nOrb2Loc,Label='xkappa_cnt') !saves the previous kappa_cnt
    call mma_Allocate(unitary_mat,nOrb2Loc,nOrb2Loc,Label='unitary_mat')
    call mma_Allocate(rotated_cmo,nBasis,nOrb2Loc,Label='rotated_cmo')
    call mma_Allocate(CMO_Ref,nBasis,nOrb2Loc,Label='CMO_Ref')

end select ! allocations


! Initialization

call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)
if (.not. Silent) write(u6,"(/A)") "MO extension before localisation:"
call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.not. Silent)


! set defaults

OldFunctional = Functional
FirstFunctional = Functional
Delta = Functional
UpMeth=" -  - "
largest=0
nDIIS=0

GEKRange = .false.
ResetGEK = .false.
mindp = 1  ! minimal number of data points for GEK construction
SORange = .true.

IterGEK = 0


! Print iteration table header.

if (.not. Silent) then
    call CWTime(C2,W2)
    TimC = C2-C1
    TimW = W2-W1
    write(u6,'(//,1X,A,/,1X,A)') '                                                                   CPU       Wall', &
    'nIter       Functional P        Delta     Gradient   Microiter   (sec)     (sec) %Screen, ndiis, largest'
end if


! ----------------------------------------------------------------------
!                           Iterations
! ----------------------------------------------------------------------

nIter = 0
Converged = .false.

do while ((nIter < nMxIter) .and. (.not. Converged))
    if (.not. Silent) call CWTime(C1,W1)

    nIter = nIter+1

    !choose between optimization methods
    select case (OptMeth)

    case (1) ! Jacobi Sweeps (2x2 rotations)

        call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
        call GetGradnorm_PM(nAtoms,nOrb2Loc,PA,GradNorm)
        call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,PctSkp)

    case (2,4,5) ! Employing NxN rotations

        call GenerateP(Ovlp,CMO,BName,nBasis,nOrb2Loc,nAtoms,nBas_per_Atom,nBas_Start,PA,Ovlp_sqrt)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Gradient(:), Hdiagvec(:)) ! gets the new gradient

        GradList(:,nIter) = -Gradient(:) ! g_i
        FuncList(nIter) = -Functional ! y_i

        ! compute standard newton raphson step
        Disp(:) = -Gradient(:)/Hdiagvec(:)

        call rescale_disp()

#       ifdef _DEBUGLISTS_
            write(u6,*) "nIter =",nIter
            call RecPrt('DispList(:,:nIter)',' ',DispList(:,:nIter),fsdim,nIter)
            call RecPrt('GradList(:,:nIter)',' ',GradList(:,:nIter),fsdim,nIter)
            call RecPrt('FuncList(:nIter)',' ',FuncList(:nIter),nIter,1)
#       endif

        if (OptMeth == 4 .or. OptMeth == 5) then ! (S)-GEK

            if (GEKRange) then
                ! still in infinitesimal limit of kappa, sampled previous point -> start GEK

                IterGEK = IterGEK + 1

                call S_GEK_localisation(nIter,IterGEK,mindp,-hdiagvec(:),fsdim,dqdq,Disp(:),UpMeth,SORange,nOrb2Loc,&
                                        usmitigation,nDIIS)

                ! undershoot mitigation
                if (usmitigation) then
                    if (Loosen%Step > One) then
                        call mma_allocate(Prev,fsdim,Label='Prev')

                        Prev(:) = DispList(:,nIter)

                        dqdq = DDot_(fsdim,Disp,1,Disp,1)*DDot_(fsdim,Prev,1,Prev,1)
                        ang = DDot_(fsdim,Prev,1,Disp,1)/sqrt(dqdq)
                        if (ang < Loosen%Thrs2) then
                            Loosen%Factor = One
                        else if (ang > Loosen%Thrs) then
                            Loosen%Factor = Loosen%Factor*Loosen%Step
                        end if

#                       ifdef _DEBUGPRINT_
                        call RecPrt('Disp',' ',Disp,fsdim,1)
                        call RecPrt('Prev',' ',Prev,fsdim,1)
                        write(u6,*) "angle(Disp,Prev) = cos^-1(",ang,")"
                        write(u6,*) "Loosen%Factor    =", Loosen%Factor
                        write(u6,*) "Loosen%Step    =", Loosen%Step
#                       endif

                        call mma_Deallocate(Prev)
                    end if ! Loosen%Step > One
                end if ! undershoot mitigation

            end if ! if in GEKRange

        end if ! NR vs GEK
        ! ---------------------------------------------------------------------------------------------------

        call rescale_disp()
        ! see if inside region fit for GEK
        call StepSizeChecks()

        ! transform disp vec to matrix
        call vec2upper_triag(kappa(:,:),nOrb2Loc,Disp(:),fsdim,.true.)
        DispList(:,nIter) = Disp(:) ! q_i

        ! update CMO
        call RotateNxN(CMO,kappa,nOrb2Loc,nBasis,kappa_cnt,xkappa_cnt,unitary_mat,rotated_CMO)
        UMatList(:,:,nIter) = unitary_mat(:,:) ! exp(-q_i) = U_i

    end select ! 2x2 or NxN rotations

#   ifdef _GETMOLDEN_
    ! choose the iteration of interest, this creates a $project.imlocal.molden file
    write (x1,fmt) nIter ! converting integer to string using a 'internal file'
    imfile = trim(NewDir)//'/imloc.'//x1//'.molden'

    ! creating files in scratch directory
    call get_intermediate_molden(CMO,nBasis,nOrb2Loc)

    ! move files from scratch dir to project dir
    call systemf("mv "//trim(WorkDir)//'/imloc '//trim(imfile),rc)
    call systemf("mv "//trim(WorkDir)//'/LocOrbIM '//trim(NewDir)//'/LocOrbIM.'//x1,rc)

#   endif


    !check if converged

    Delta = Functional-OldFunctional
    OldFunctional = Functional
    if (.not. Silent) then
        call CWTime(C2,W2)
        TimC = C2-C1
        TimW = W2-W1
        write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),3X,A6,1X,2(1X,F9.1),1X,F7.2,1X,I5,1X,ES12.4)') &
            nIter,Functional,Delta,GradNorm,UpMeth,TimC,TimW,PctSkp,nDIIS,largest
    end if

    Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)

end do !Iterations


#ifdef _DEBUGPRINT_
write(u6,'(/A)') 'Check the orthonormality of the orbitals'
write(u6,*) '========================================'
call dgemm_('T','N',nOrb2Loc, nBasis, nBasis,One, CMO, nBasis,Ovlp, nBasis,Zero, CtS, nOrb2Loc)
call dgemm_('N','N',nOrb2Loc, nOrb2Loc, nBasis,One,CtS, nOrb2Loc,CMO, nBasis,Zero,CtSC, nOrb2Loc)
call RecPrt("C^T*S*C =",' ',CtSC,nOrb2Loc, nOrb2Loc)
#endif

! Print convergence message.

if (.not. Silent) then

    write(u6,"(/A)") "MO extension after localisation:"
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)

    if (.not. Converged) then
        write(u6,'(/,A,I4,A)') 'No convergence after',nIter,' iterations.'
    else
        write(u6,'(/,A,I4,A)') 'Convergence after',nIter,' iterations.'
        write(u6,*)
        write(u6,'(A,I8)') 'Number of localised orbitals  : ',nOrb2loc
        write(u6,'(A,ES20.10)') 'Value of P before localisation: ',FirstFunctional
        write(u6,'(A,ES20.10)') 'Value of P after localisation : ',Functional
        write(u6,*) "Localisation converged in ",nIter
        write(u6,*) "Last deltaP",Delta
        write(u6,*) "Last GradNorm",GradNorm
    end if
end if

!call Prpt()


! deallocations

select case(OptMeth)
case(1)
    call mma_Deallocate(PACol)
case(2,3,4,5)
    call mma_Deallocate(Gradient)
    call mma_Deallocate(Hdiag)
    call mma_Deallocate(kappa)

    call mma_Deallocate(kappa_cnt)
    call mma_Deallocate(xkappa_cnt)
    call mma_Deallocate(unitary_mat)
    call mma_Deallocate(rotated_CMO)
    call mma_Deallocate(CMO_Ref)

    call mma_Deallocate(FuncList)
    call mma_Deallocate(UmatList)
    call mma_Deallocate(GradList)
    call mma_Deallocate(DispList)
    call mma_Deallocate(disp)
    call mma_Deallocate(Hdiagvec)
end select

call mma_Deallocate(Ovlp_sqrt)

contains
subroutine rescale_disp()

DD=Sqrt(dot_product(Disp(:),Disp(:)))
!Thr= 0.5E0_wp * Pi
Thr= 0.25E0_wp * Pi
If (DD>=Thr)Then
#   ifdef _DEBUGPRINT_
    Write(u6,*) 'Rescale Kappa(:)'
#   endif
Disp(:) = (Thr/DD)*Disp(:)
End If

end subroutine rescale_disp

subroutine StepSizeChecks()
        ! if previous step suggestion was out of GEKRange
        if (ResetGEK) then
            UpMeth=" -  - "
            IterGEK = 0
            nDIIS = 0
            ResetGEK = .false.
        end if

        ! check if matrix elements are > 0.01
        large_elements = 0
        do i=1,fsdim
            if (abs(Disp(i)) > 0.005_wp) then
                large_elements = large_elements + 1
            else
                large_elements = large_elements
            end if
        end do
        maxel = maxloc(abs(Disp),1)
        largest = Disp(maxel)

        ! all elements of kappa are small enough to use this disp as coordinate for building the GEK model
        if (large_elements == 0) then
            GEKRange = .true.
        else if (large_elements /= 0 .and. GEKRange .and. IterGEK > 0) then
            ! leave GEK and go back to NR if steps are too large
            !write(u6,*) "resetting GEK due to large step:",largest

            !Disp(:) = 0.1*Disp(:)
            ResetGEK = .true.
            GEKRange = .false.
        else
            GEKRange = .false.
        end if

#       ifdef _DEBUGPRINT_
        write(u6,*) "kappa elements > 0.01 =",large_elements
        write(u6,*) "largest element =", Disp(maxel)
#       endif
end subroutine StepSizeChecks

end subroutine PipekMezey_Iter
