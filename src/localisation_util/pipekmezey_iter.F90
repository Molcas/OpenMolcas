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
!               2026, Lila Zapp                                        *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _FORCEGEKRANGE_
!#define _TESTNUMERICALLY_
!#define _JSINSP_
subroutine PipekMezey_Iter(Functional,CMO,PA,nBasis,nOrb2Loc,Converged)
! Author: T.B. Pedersen
!
! Based on the original routines by Y. Carissan.
! 2026, Lila Zapp (opt methods & Lowdin framework)

use Index_Functions, only: nTri_Elem
use Localisation_globals, only: AnalyseLoc, bias, BName, ChargeType, Debug, Debug, DispList, FuncList, GEKThr_Grad, GEKThr_Kappa, &
                                getIMmldn, GradList, inpOptMeth, kappa_cnt, MoldMod, nAtoms, nBas_per_Atom, nBas_Start, nMxIter, &
                                OptMeth, Ovlp, Ovlp_sqrt, posel, ThrGrad, Thrs, xkappa_cnt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc
real(kind=wp), intent(out) :: Functional, PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
logical(kind=iwp), intent(out) :: Converged
integer(kind=iwp) :: fsdim, IterGEK, large_elements, lSCR, maxel, nDIIS, nIter, npos
#ifdef _JSINSP_
integer(kind=iwp) :: kl
#endif
real(kind=wp) :: a, alpha, b, best_eta, C1, C2, DD, Delta, dqdq, eta1, eta2, FirstFunctional, GradNorm, largest, OldFunctional, &
                 P_eta0, P_eta1, P_eta2, PctSkp, scalingfac, StepNorm, Thr, TimC, TimW, W1, W2
logical(kind=iwp) :: GEKRange, linesearch, ResetGEK, switched
character(len=6) :: UpMeth
logical(kind=iwp), parameter :: debug_lowd = .false., modHess = .true., trafoPA = .false.
real(kind=wp), allocatable :: CMO_Ref(:,:), Disp(:), Gradient(:), hdiagvec(:), kappa(:,:), Ovlp_aux(:,:), PACol(:,:), &
                              rotated_CMO(:,:), SCR(:), SearchDir(:), Umat(:,:)
#ifdef _TESTNUMERICALLY_
real(kind=wp), allocatable :: NumGrad(:)
#endif
#if defined(_TESTNUMERICALLY_) || defined(_JSINSP_)
real(kind=wp), allocatable :: Hessian(:)
#endif

call CWTime(C1,W1)

! remember input choice of optimization method
OptMeth = inpOptMeth

! print info for GEK settings
if ((OptMeth == 4) .or. (OptMeth == 5) .or. (OptMeth == 6)) then
  write(u6,*)
  write(u6,*) 'Settings for (S)-GEK optimization'
  write(u6,'(1X,A,ES12.4)') 'GEKThrKappa          :',GEKThr_Kappa
  write(u6,'(1X,A,ES12.4)') 'GEKThrGrad           :',GEKThr_Grad
  write(u6,'(1X,A,ES12.4)') 'bias                 :',bias
end if

! ensure that the initial orbitals are orthonormal
call OrthoCheck(CMO,nOrb2Loc,nBasis)

! to allow property printing later
! aim: compute the locality of the MOs using second and fourth moment orbital spread (not completed)
!call Put_cArray('Relax Method','LOCALIS ',8)

! if the Lowdin charge framework is requested instead of Mulliken
call mma_allocate(Ovlp_sqrt,nBasis,nBasis,Label='S^{1/2}')
if (ChargeType == 2) then
  call mma_allocate(Ovlp_aux,nBasis,nBasis,Label='S^{-1/2}')
  lSCR = 2*nBasis**2+nTri_Elem(nBasis)

  if (debug_lowd) call RecPrt('S before taking the sqrt',' ',Ovlp,nBasis,nBasis)

  call mma_allocate(SCR,lSCR,Label='SCR')
  call SQRTMT(Ovlp,nBasis,1,Ovlp_sqrt,Ovlp_aux,SCR)
  call mma_Deallocate(SCR)

  if (debug_lowd) then
    call RecPrt('S^{1/2}',' ',Ovlp_sqrt,nBasis,nBasis)
    call RecPrt('S after taking the sqrt',' ',Ovlp,nBasis,nBasis)
    Ovlp_aux(:,:) = Zero ! I want to reuse it
    call dgemm_('N','N',nBasis,nBasis,nBasis,One,Ovlp_sqrt,nBasis,Ovlp_sqrt,nBasis,Zero,Ovlp_aux,nBasis)
    call RecPrt('S^{1/2}*S^{1/2}',' ',Ovlp_aux,nBasis,nBasis) ! should be same as S
  end if

  call mma_deallocate(Ovlp_aux)
end if

! number of orbital pairs to localise / number of parameters (fullspace dimension)
fsdim = nTri_Elem(nOrb2Loc-1)

! allocations
! -----------
call mma_Allocate(Hdiagvec,fsdim,Label='Hdiagvec')

!select case (InpOptMeth)

!  case (1)

!call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')

!  case (2,3,4,5,6)

! allocating matrices for NxN optimizations

call mma_Allocate(kappa,nOrb2Loc,nOrb2Loc,Label='kappa')
call mma_Allocate(Gradient,fsdim,Label='Gradient')
call mma_Allocate(posel,fsdim,Label='posel')
call mma_Allocate(DispList,fsdim,nMxIter,Label='DispList')  ! kappa matrices
call mma_allocate(Disp,fsdim,Label='Disp')
call mma_allocate(SearchDir,fsdim,Label='SearchDir')
call mma_Allocate(GradList,fsdim,nMxIter,Label='GradList')
call mma_Allocate(FuncList,nMxIter,Label='FuncList')
call mma_Allocate(kappa_cnt,nOrb2Loc,nOrb2Loc,Label='kappa_cnt') != kappa^cnt
call mma_Allocate(xkappa_cnt,nOrb2Loc,nOrb2Loc,Label='xkappa_cnt') !saves the previous kappa_cnt
call mma_Allocate(Umat,nOrb2Loc,nOrb2Loc,Label='Umat')
call mma_Allocate(rotated_cmo,nBasis,nOrb2Loc,Label='rotated_cmo')
call mma_Allocate(CMO_Ref,nBasis,nOrb2Loc,Label='CMO_Ref')

Umat(:,:) = Zero
FuncList(:) = Zero
Kappa(:,:) = Zero
Gradient(:) = Zero
GradList(:,:) = Zero
DispList(:,:) = Zero
posel(:) = 0
if ((OptMeth == 1) .or. (OptMeth == 6)) call mma_Allocate(PACol,nOrb2Loc,2,Label='PACol')

!  case default
!    write(u6,*) 'ERROR: The chosen opt method is not implemented for localisation'
!    call Abend()

!end select ! allocations

! Initialization

call GenerateP(CMO,nBasis,nOrb2Loc,nAtoms,PA)

select case (AnalyseLoc)
  case (0,1)
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
  case (2)
    write(u6,'(/A)') 'Orbital extension d_i before localisation:'
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)
  case default
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
end select

!select case (OptMeth)

!  case (1,6)
call GetGradnorm_PM(nAtoms,nOrb2Loc,PA,GradNorm)

!  case (2,3,4,5)
call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Gradient) ! gets the initial gradient
FuncList(1) = -Functional
GradList(:,1) = -Gradient(:)

#ifdef _TESTNUMERICALLY_
call mma_allocate(NumGrad,fsdim,Label='NumGrad')
call GetNumGrad_PM(CMO,nOrb2Loc,nBasis,fsdim,NumGrad,.true.)
call mma_deallocate(NumGrad)

call mma_allocate(Hessian,fsdim,fsdim,Label='Hessian')
call GetHess_PM(nAtoms,nOrb2Loc,PA,fsdim,Hessian,CMO,nBasis)
call RecPrt('full analytical Hessian','(21F10.6)',Hessian,fsdim,fsdim)
call mma_deallocate(Hessian)
#endif

!  case default
!    write(u6,*) 'ERROR: for the selected OptMeth, there exists no initialisation'
!    call Abend()

!end select

FirstFunctional = Functional
Delta = Functional
largest = Zero
nDIIS = 0
OldFunctional = Zero
scalingfac = One
GEKRange = .false.
ResetGEK = .false.
switched = .false.
IterGEK = -1
PctSkp = 0
Converged = .false.

! Print iteration table header.
! -----------------------------

write(u6,*)
write(u6,*) 'npos = number of positive Hessian diagonal elements'

call CWTime(C2,W2)
TimC = C2-C1
TimW = W2-W1
nIter = 0
! initial information (Iteration = 0)
select case (InpOptMeth)
  case (1)
    UpMeth = 'JS  - '
    write(u6,'(//,1X,A,/,1X,A)') '                                                                 CPU       Wall', &
                                 'nIter       Functional P        Delta     Gradient   Method     (sec)     (sec)   npos  %Screen'
  case (2,4,5)
    UpMeth = 'NR   0'
    write(u6,'(//,1X,A,/,1X,A)') '                                                                 CPU       Wall', &
                                 'nIter       Functional P        Delta     Gradient   Method     (sec)     (sec)   npos  dispnorm'
  case (3)
    UpMeth = 'GA  - '
    write(u6,'(//,1X,A,/,1X,A)') '                                                                 CPU       Wall', &
                                 'nIter       Functional P        Delta     Gradient   Method     (sec)     (sec)  npos  dispnorm'
  case (6)
    UpMeth = 'JS  - '
    write(u6,'(//,1X,A,/,1X,A)') '                                                                 CPU       Wall', &
                                 'nIter       Functional P        Delta     Gradient   Method     (sec)     (sec)  npos  '// &
                                 '%Screen/dispnorm'
  case default
    write(u6,*) 'ERROR: The chosen opt method is not implemented for localisation'
    call Abend()
end select

! ----------------------------------------------------------------------
!                           Iterations
! ----------------------------------------------------------------------

linesearch = .false.

do while ((nIter < nMxIter) .and. (.not. Converged))
  call CWTime(C1,W1)

  nIter = nIter+1
  if (Debug) write(u6,*) 'nIter = ',nIter

  call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
  call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Gradient)
  call GetHdiag_PM(nAtoms,nOrb2Loc,PA,Hdiagvec,npos,gradnorm,modHess)

  if (inpOptMeth == 6) then
    if ((npos /= 0) .and. (.not. switched)) then
      ! request to start with one Jacobi Sweep, then switch to NR (6) or GEK (7)
      OptMeth = 1
    else
      OptMeth = 5
      switched = .true.
    end if
  end if

  ! choose between optimization methods
  select case (OptMeth)

    case (1) ! Jacobi Sweeps (2x2 rotations)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! PAIR ROTATIONS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      UpMeth = 'JS  -'

      call RotateOrb(CMO,PACol,nBasis,nAtoms,PA,nOrb2Loc,BName,nBas_per_Atom,nBas_Start,PctSkp)

    case (2,3,4,5,6)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! N X N ROTATIONS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! EVALUATE FUNCTIONAL, GRADIENT, HESSIAN DIAGONAL AT CURRENT POINT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FuncList(nIter) = -Functional ! y_i
      GradList(:,nIter) = -Gradient(:) ! g_i

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! GRADIENT ASCENT STEP
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (OptMeth == 3) then
        ! Gradient Ascent with naive line search
        !Disp(:) = Gradient(:)/gradnorm

        linesearch = .true.
        !SearchDir(:) = Gradient(:)/gradnorm
        SearchDir(:) = Gradient(:)

        if (linesearch) then
          alpha = One
          call naiveLineSearch(SearchDir,best_eta,alpha)
          write(u6,*) 'best_eta',best_eta
          UpMeth = 'GA +LS'
        else
          best_eta = One
          UpMeth = 'GA  -'
        end if

        Disp(:) = best_eta*SearchDir(:)

      else
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! NEWTON RAPHSON STEP
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SearchDir(:) = -Gradient(:)/Hdiagvec(:)

#       ifdef _JSINSP_
        ! experimental code aiming to escape stationary points while the hessian still has positive eigenvalues
        call mma_allocate(Hessian,fsdim,fsdim,Label='Hessian')
        call GetHess_PM(nAtoms,nOrb2Loc,PA,fsdim,Hessian,CMO,nBasis)
        !call RecPrt('full analytical Hessian','(6F10.6)',Hessian,fsdim,fsdim)
        call mma_deallocate(Hessian)

        !call RecPrt('SearchDir',' ',SearchDir,fsdim,1)
        !call RecPrt('Hdiagvec',' ',Hdiagvec,fsdim,1)
        !call RecPrt('Gradient',' ',Gradient,fsdim,1)

        if ((npos > 0) .and. (gradnorm > 1.0e-2_wp)) then
          do kl=1,fsdim
            if ((posel(kl) == 1) .and. (abs(Gradient(kl)) < 1.0e-2_wp)) then
              !SearchDir(kl) = Pi/Hdiagvec(kl)
              SearchDir(kl) = Pi
              write(u6,*) 'kl,posel(kl),Gradient(kl),Hessian(kl),SearchDir(kl)',kl,posel(kl),Gradient(kl),Hdiagvec(kl),SearchDir(kl)
            end if
          end do
        end if
        !call RecPrt('SearchDir',' ',SearchDir,fsdim,1)
#       endif

        if (linesearch) then
          alpha = One
          call naiveLineSearch(SearchDir,best_eta,alpha)
          UpMeth = 'NR +LS'
        else
          best_eta = One
          UpMeth = 'NR   0'
        end if

        Disp(:) = best_eta*SearchDir(:)
      end if

      ! keep norm of kappa matrix below pi/4
      call rescale_disp(Disp)

      ! see if inside region fit for GEK
      if (OptMeth > 3) call StepSizeChecks()

#     ifdef _FORCEGEKRANGE_
      call force_GEKRange()
#     endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! IF IN GEK RANGE: Build subspace and get back optimized Disp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((OptMeth == 4) .or. (OptMeth == 5) .or. (OptMeth == 6)) then ! (S)-GEK
        if (GEKRange) then

          IterGEK = IterGEK+1

          call S_GEK_localisation(nIter,IterGEK,fsdim,dqdq,Disp,UpMeth,nOrb2Loc,nDIIS,-hdiagvec,CMO,nBasis,PA,nAtoms)

        end if ! if in GEKRange

      end if ! NR vs GEK
      ! ---------------------------------------------------------------------------------------------------
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! TAKE THE GEK STEP; REPLACE THE NR DATA IN THE LISTS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! keep norm of kappa matrix below pi/4
      call rescale_disp(Disp)

      ! see if inside region fit for GEK
      call StepSizeChecks()

      ! transform disp vec to matrix
      call vec2upper_triag(kappa,nOrb2Loc,Disp,fsdim,.true.)

      ! get U=exp(-kappa) and U_inv=exp(kappa)
      call expkap_localisation(kappa,nOrb2Loc,Umat)
      ! rotate MOs as rotated_CMO = CMO * exp(-kappa)
      call RotateNxN(CMO,nOrb2Loc,nBasis,Umat,rotated_CMO)
      ! update <s|PA|t>

      if (trafoPA) then
        call transformPA(PA,nOrb2Loc,Umat,.true.)
      else
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
      end if

#     ifdef _DEBUGPRINT_
      call RecPrt('Gradient',' ',Gradient,fsdim,1)
      call RecPrt('Hdiagvec',' ',Hdiagvec,fsdim,1)
      write(u6,*) 'After GEK procedure and step scaling'
      call RecPrt('Disp','(F18.9)',Disp,fsdim,1)
      call RecPrt('Unitary Mat','(10F18.9)',Umat,norb2loc,norb2loc)
#     endif

      DispList(:,nIter) = Disp(:) ! q_i
      CMO(:,:) = rotated_CMO(:,:)
      Stepnorm = sqrt(dot_product(Disp,Disp))
  end select ! 2x2 or NxN rotations

  if (getIMmldn .and. (mod(nIter,MoldMod) == 0)) call get_intermediate_molden(nIter)

  !check if converged
  Delta = Functional-OldFunctional
  OldFunctional = Functional
  call CWTime(C2,W2)
  TimC = C2-C1
  TimW = W2-W1
  !write(u6,*)
  select case (OptMeth)
    case (1)
      write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),3X,A6,1X,2(F9.3,1X),I5,1X,F8.2)') &
        nIter-1,Functional,Delta,GradNorm,UpMeth,TimC,TimW,npos,PctSkp
    case (2,3,4,5,6)
      write(u6,'(1X,I5,1X,F18.8,2(1X,ES12.4),3X,A6,1X,2(F9.3,1X),I5,1X,ES12.4)') &
        nIter-1,Functional,Delta,GradNorm,UpMeth,TimC,TimW,npos,StepNorm
    case default
      write(u6,*) 'ERROR: The chosen opt method is not implemented for localisation'
      call Abend()
  end select

  Converged = (GradNorm <= ThrGrad) .and. (abs(Delta) <= Thrs)

end do !Iterations

call OrthoCheck(CMO,nOrb2Loc,nBasis)

#ifdef _TESTNUMERICALLY_
call mma_allocate(NumGrad,fsdim,Label='NumGrad')
call GetNumGrad_PM(CMO,nOrb2Loc,nBasis,fsdim,NumGrad,.true.)
call mma_deallocate(NumGrad)

call mma_allocate(Hessian,fsdim,fsdim,Label='Hessian')
call GetHess_PM(nAtoms,nOrb2Loc,PA,fsdim,Hessian,CMO,nBasis)
call RecPrt('full analytical Hessian','(6F10.6)',Hessian,fsdim,fsdim)
call mma_deallocate(Hessian)
#endif

! Print convergence message.
! --------------------------
call GenerateP(CMO,nBasis,nOrb2Loc,nAtoms,PA)
select case (AnalyseLoc)
  case (0)
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)
  case (1,2)
    write(u6,'(/A)') 'Orbital extension d_i after localisation:'
    call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.true.)
end select

if (.not. Converged) then
  write(u6,'(/,A,I4,A)') 'No convergence after',nIter-1,' iterations.'
else
  write(u6,'(/,A,I4,A)') 'Convergence after',nIter-1,' iterations.'
  write(u6,*)
  write(u6,'(A,I8)') 'Number of localised orbitals  : ',nOrb2loc
  write(u6,'(A,ES20.10)') 'Value of P before localisation: ',FirstFunctional
  write(u6,'(A,ES20.10)') 'Value of P after localisation : ',Functional
end if

!call Prpt()

call Add_Info('LOC_ITER',[real(nIter,kind=wp)],1,8)

! deallocations
! -------------
call mma_Deallocate(PACol,safe='*')
call mma_Deallocate(Hdiagvec)
!select case (InpOptMeth)
!  case (2,3,4,5,6)
call mma_Deallocate(Gradient)
call mma_Deallocate(kappa)
call mma_Deallocate(posel)

call mma_Deallocate(kappa_cnt)
call mma_Deallocate(xkappa_cnt)
call mma_Deallocate(Umat)
call mma_Deallocate(rotated_CMO)
call mma_Deallocate(CMO_Ref)

call mma_Deallocate(FuncList)
call mma_Deallocate(GradList)
call mma_Deallocate(DispList)
call mma_Deallocate(Disp)
call mma_Deallocate(SearchDir)

!end select

call mma_Deallocate(Ovlp_sqrt)

contains

subroutine rescale_disp(Disp)
  ! the current exp(-kappa) routine cannot handle large displacements,
  ! so we rescale the norm of kappa down to 0.99 Pi

  real(kind=wp), intent(inout) :: Disp(fsdim)

  DD = sqrt(dot_product(Disp,Disp))
  Thr = 0.99_wp*Pi
  if (DD >= Thr) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Rescale Kappa(:)'
#   endif
    scalingfac = Thr/DD
    Disp(:) = scalingfac*Disp(:)
  else
    scalingfac = One
  end if

end subroutine rescale_disp

#ifdef _FORCEGEKRANGE_
subroutine force_GEKRange()
  ! scales stepsize down to fulfill gekthr_kappa criterium; ignores the grad criterium by overwriting it

  gekthr_grad = 100.0_wp
  Thr = 0.99_wp*gekthr_kappa
  maxel = maxloc(abs(Disp),1)
  largest = Disp(maxel)
  if (abs(largest) > Thr) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Rescale Disp(:)'
    write(u6,*) 'largest',largest
    write(u6,*) 'Disp(:) =',Disp(:),'Thr/abs(largest)*Disp(:) = ',Thr/abs(largest)*Disp(:)
#   endif
    Disp(:) = Thr/abs(largest)*Disp(:)
  end if

end subroutine force_GEKRange
#endif

subroutine StepSizeChecks()
  ! this subroutine checks if the Range is reached, in which data points can be collected for the GEK
  ! based on a) GEKThr_Kappa, avoiding that single elements of the displacement vector kappa are too large
  ! which is required because the coordinate mapping assumes linear coordinate behavior; and
  ! b) GEKThr_Grad, which is an additional criterion (should be the preferred one as soon as the coordinate
  ! mapping is done exactly and not via linear approximation);
  ! both criteria have to be fulfilled for the GEK to start/continue sampling
  ! if one of them is not fulfilled, the GEK will be resetted and switched back to NR

  integer(kind=iwp) :: i

  ! if previous step suggestion was out of GEKRange
  if (ResetGEK .and. (nDIIS /= 0)) then
    !write(u6,*) 'Resetting GEK'
    IterGEK = 0
    nDIIS = 0
    ResetGEK = .false.
    OptMeth = InpOptMeth
  end if

  ! check if matrix elements are > 0.01
  large_elements = 0
  do i=1,fsdim
    if (abs(Disp(i)) > gekthr_kappa) then
      large_elements = large_elements+1
    else
      large_elements = large_elements
    end if
  end do
  maxel = maxloc(abs(Disp),1)
  largest = Disp(maxel)

  ! all elements of kappa are small enough to use this disp as coordinate for building the GEK model
  if ((large_elements == 0) .and. (GradNorm < gekthr_grad) .and. (npos == 0)) GEKRange = .true.

  if (GEKRange) then
    if (Gradnorm >= gekthr_grad) then
      if (Debug) write(u6,'(A,ES18.8)') 'Reset GEK in next iteration because Grad >',gekthr_grad
      ResetGEK = .true.
      GEKRange = .false.
    end if
    if (large_elements /= 0) then
      if (Debug) &
        write(u6,'(A,ES18.8)') 'Reset GEK in next iteration because largest element of the kappa matrix is above',gekthr_kappa
      ResetGEK = .true.
      GEKRange = .false.
    end if
  end if

# ifdef _DEBUGPRINT_
  write(u6,*) 'kappa elements > 0.01 =',large_elements
  write(u6,*) 'largest element =',Disp(maxel)
# endif

end subroutine StepSizeChecks

subroutine P_of_eta(Disp,Functional)

  real(kind=wp), intent(in) :: Disp(fsdim)
  real(kind=wp), intent(out) :: Functional

  call vec2upper_triag(kappa,nOrb2Loc,Disp,fsdim,.true.)
  call expkap_localisation(kappa,nOrb2Loc,Umat)
  call RotateNxN(CMO,nOrb2Loc,nBasis,Umat,rotated_CMO)
  call transformPA(PA,nOrb2Loc,Umat,.true.)
  call ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,.false.)

end subroutine P_of_eta

subroutine naiveLineSearch(SearchDir,best_eta,alpha)

  real(kind=wp), intent(inout) :: SearchDir(fsdim), alpha
  real(kind=wp), intent(out) :: best_eta

  call rescale_disp(SearchDir)

  a = -One

  do while (a < Zero)
    write(u6,*) 'alpha = ',alpha
    eta1 = Half*alpha
    eta2 = alpha

    P_eta0 = Functional
    call P_of_eta(Zero*SearchDir(:),P_eta0) ! just checking, this equals current functional
    call P_of_eta(eta1*SearchDir(:),P_eta1)
    call P_of_eta(eta2*SearchDir(:),P_eta2)

    write(u6,*) 'P_eta0 =',P_eta0
    write(u6,*) 'P_eta1 =',P_eta1
    write(u6,*) 'P_eta2 =',P_eta2
    !write(u6,*) 'P_eta0,P_eta1, P_eta2 =',P_eta0,P_eta1,P_eta2
    b = ((P_eta2-P_eta0)*eta1**2-(P_eta1-P_eta0)*eta2**2)/(eta1*eta2*(eta1-eta2))
    a = (P_eta2-P_eta0)/(eta2**2)-b/eta2
    write(u6,*) 'a =',a,'b=',b,'c=',P_eta0
    ! we want the maximum, not the minimum, so f'(eta)=2a must be positive

    best_eta = -b/(Two*a)
    !write(u6,*) 'best_eta =',best_eta
    if (a < Zero) then
      alpha = Two*alpha
      write(u6,*) 'LS lead to minimum, rescaling kappa now'
    end if
  end do

end subroutine naiveLineSearch

end subroutine PipekMezey_Iter
