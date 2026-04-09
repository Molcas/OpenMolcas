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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************
!
! This module is used when we need to consider rotations between
! internal states, which are usually state-averaged RASSCF states. The
! standard SA-RASSCF (CASSCF) is invariant with respect to rotations
! within internal states, but not for the following example:
!   1. unequally weighted SA-RASSCF
!   2. SA-RASSCF/PCM (unless state-averaging and density that polarizes
!      ASCs are the same (def_solv = 3))
!   3. CASPT2/RASPT2
! In some cases, the Lagrangian multipliers relevant to the rotations
! can be obtained without iteration. For instance, for 3, if the
! reference SA-CASSCF is invariant, the multipliers can be immediately
! obtained. Or, for 2, if the density that polarizes ASCs are equally
! state-averaged, the same applies. For other cases, they should be
! obtained iteratively.
!
! (InvSCF, InvEne) = (.true., .true.)
!   -> conventional SA-CASSCF/RASSCF and SA-CASSCF/PCM (with def_solv = 3)
!      There is no need to consider the internal state rotations
! (InvSCF, InvEne) = (.true., .false.)
!   -> CASPT2 (w/o PCM), SA-CASSCF/PCM (with def_solv = 4)
!      Internal state rotations are non-iteratively evaluated by rotating
!      the 1e and 2e density matrix in rhs_sa and rhs_nac through the
!      SLag matrix.
! (InvSCF, InvEne) = (.false., .false.)
!   -> anything based on SA-CASSCF/PCM (with def_solv /= 3,4)
!      Some internal rotations may be evaluated non-iteratively as in
!      (InvSCF, InvEne) = (.true., .false.). This comes from the
!      rotation in MS-CASPT2 and should not be considered as rotation
!      parameters to be optimized, since it is not iteratively
!      determined during energy.
!      The other internal rotations come from the overlap between the
!      original CI and CI derivative. This is optimized during Z-vector.
! (InvSCF, InvEne) = (.false., .true.)
!   -> unequally (including dynamically) weighted SA-MCSCF w/o PCM
!      The initial residue of the internal rotation is zero, though.

module ISRotation

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none

! Whether the "target" energy, specified by RLXROOT, is invariant or
! non-invariant wrt internal state rotations. The conventional
! SA-MCSCF is invariant, so InvEne = .true. For the above cases, the
! energy is non-invariant, so InvEne = .false. (for instance, CASPT2
! or SA-MCSCF with PCM etc.). This flag only controls the evaluation
! of the initial internal rotation parameter.
logical(kind=iwp) :: InvEne = .true.
! Whether the reference "state-averaged" SCF energy is invariant or
! not with respect to rotations between internal states. If they are
! invariant (InvSCF = .true.), the rotation parameters can be
! non-iteratively determined. If non-invariant (InvSCF = .false.),
! they have to be iteratively determined during the Z-vector
! procedure. The rotation parameters for InvSCF = .true. can be
! iteratively determined as well; this will end up with the same
! result with InvSCF = .false. See the above note for details.
logical(kind=iwp) :: InvSCF = .true.
! Whether we consider the off-diagonal block (O-S and C-S) couplings.
! Purely development purpose
logical(kind=iwp) :: IntRotOff = .true.
! unequal state-average?
logical(kind=iwp) :: unequal_SA = .false.
! InvSol is true only if the ASCs are polarized by the state-averaged
! density matrix. If .not.InvSCF, InvSol is always false, but the
! converse may not be true.
logical(kind=iwp) :: InvSol = .true.
! Whether state rotations are scaled with the weight factor.
! This option should not affect the computed results but I'm not sure.
! Purely development purpose
logical(kind=iwp) :: ScalWeight = .false.

type ISR_param
  real(kind=wp), allocatable :: Ap(:,:)   ! results of A*p
  real(kind=wp), allocatable :: p(:,:)    ! trial vector during CG
  real(kind=wp), allocatable :: Rvec(:,:) ! initial residue (RHS)
  real(kind=wp), allocatable :: prec(:,:) ! preconditioned something
  real(kind=wp), allocatable :: Xvec(:,:) ! solution
  real(kind=wp), allocatable :: Pvec(:,:) ! P vector (for CGS)
  real(kind=wp), allocatable :: Qvec(:,:) ! Q vector (for CGS)
  real(kind=wp), allocatable :: Uvec(:,:) ! U vector (for CGS)
  real(kind=wp), allocatable :: R0(:,:)   ! initial residue (for CGS)
end type ISR_param

type(ISR_param) :: ISR

logical(kind=iwp), pointer :: do_RF

contains

!-----------------------------------------------------------------------

subroutine ISR_init(iPL,do_RF_,def_solv)

  !use DWSol, only: DWSCF
  use cgs_mod, only: CGS
  use input_mclr, only: nRoots, PT2, weight

  integer(kind=iwp), intent(in) :: iPL, def_solv
  logical(kind=iwp), intent(in), target :: do_RF_
  integer(kind=iwp) :: iR

  !! Check whether we need to consider internal state rotations
  !! iteratively or non-interatively
  InvEne = .true.
  unequal_SA = .false.

  if (.not. InvSCF) InvEne = .false.
  if (PT2 .and. (nRoots > 0)) InvEne = .false.
  if (do_RF_ .and. (def_solv /= 3)) InvEne = .false.
  do iR=2,nRoots
    if (weight(1) /= weight(iR)) then
      unequal_SA = .true.
      InvSCF = .false.
      if (do_RF_ .and. (def_solv == 3)) InvEne = .false.
    end if
  end do

  do_RF => do_RF_

  !! If InvEnv = .false., we need internal state rotations
  !! Rvec may be allocated for all cases
  call mma_allocate(ISR%Rvec,nRoots,nRoots,Label='ISR%Rvec')
  ISR%Rvec(:,:) = Zero
  if (.not. InvSCF) then
    call mma_allocate(ISR%Ap,nRoots,nRoots,Label='ISR%Ap')
    call mma_allocate(ISR%p,nRoots,nRoots,Label='ISR%p')
    call mma_allocate(ISR%prec,nRoots,nRoots,Label='ISR%prec')
    call mma_allocate(ISR%Xvec,nRoots,nRoots,Label='ISR%Xvec')

    ISR%Ap(:,:) = Zero
    ISR%p(:,:) = Zero
    ISR%prec(:,:) = Zero
    ISR%Xvec(:,:) = Zero
  end if

  !! The rest is for CGS
  !! actually, if InvSCF = .false., it is usually better to use CGS
  CGS = CGS .or. do_RF ! .or. DWSCF%do_DW
  if (CGS) then

    call mma_allocate(ISR%Pvec,nRoots,nRoots,Label='ISR%Pvec')
    call mma_allocate(ISR%Qvec,nRoots,nRoots,Label='ISR%Qvec')
    call mma_allocate(ISR%Uvec,nRoots,nRoots,Label='ISR%Uvec')
    call mma_allocate(ISR%R0,nRoots,nRoots,Label='ISR%R0')

    ISR%Pvec(:,:) = Zero
    ISR%Qvec(:,:) = Zero
    ISR%Uvec(:,:) = Zero
    ISR%R0(:,:) = Zero
  end if

  if (iPL >= 2) then
    if (.not. InvEne) write(u6,'(1X,A)') 'Target energy is non-invariant with respect to internal state rotations'
    if (.not. InvSCF) write(u6,'(1X,A)') 'Internal rotations are iteratively determined'
    if (.not. InvSol) write(u6,'(1X,A)') 'Solvent ESP is non-invariant with respect to internal state rotations'
    if (CGS) write(u6,'(1X,A)') 'Preconditioned conjugate gradient squared (CGS) method is used'
    if (.not. InvEne .or. (.not. InvSCF) .or. (.not. InvSol) .or. CGS) write(u6,*)
  end if

end subroutine ISR_init

!-----------------------------------------------------------------------

subroutine ISR_final()

  use cgs_mod, only: CGS

  call mma_deallocate(ISR%Rvec)

  !! Not quite safe for unequally weighted SA-RASSCF...
  if (.not. InvSCF) then
    call mma_deallocate(ISR%Ap)
    !call mma_deallocate(ISR%p) ! deallocated in ISR_final2
    call mma_deallocate(ISR%prec)
    call mma_deallocate(ISR%Xvec)
  end if

  if (CGS) then
    call mma_deallocate(ISR%Pvec)
    call mma_deallocate(ISR%Qvec)
    call mma_deallocate(ISR%Uvec)
    call mma_deallocate(ISR%R0)
  end if

end subroutine ISR_final

!-----------------------------------------------------------------------

subroutine ISR_final2()

  if (.not. InvSCF) call mma_deallocate(ISR%p)

end subroutine ISR_final2

!-----------------------------------------------------------------------

subroutine ISR_RHS(CI,CIDER)

  use input_mclr, only: ERASSCF, ncsf, nRoots, State_Sym

  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots), CIDER(ncsf(State_Sym),nRoots)
  integer(kind=iwp) :: i, j
  real(kind=wp) :: scal
  real(kind=wp), external :: ddot_

  ! Construct the RHS (with minus) of state rotation

  do i=1,nRoots
    do j=1,i-1
      scal = DDot_(ncsf(State_Sym),CI(:,i),1,CIDER(:,j),1)-DDot_(ncsf(State_Sym),CI(:,j),1,CIDER(:,i),1)
      if (InvSCF) then
        ! non-iterative internal state rotations, if invariant
        ISR%Rvec(i,j) = ISR%Rvec(i,j)+scal/(ERASSCF(j)-ERASSCF(i))
      else
        ! initial residue (for iterative solution)
        ISR%Rvec(i,j) = ISR%Rvec(i,j)+scal
      end if
    end do
  end do

  !write(6,*) 'initial state rotation'
  !call sqprt(isr%rvec,nroots)

end subroutine ISR_RHS

!-----------------------------------------------------------------------

subroutine ISR_projection(CI,CIDER)

  use input_mclr, only: ncsf, nRoots, State_Sym

  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots)
  real(kind=wp), intent(inout) :: CIDER(ncsf(State_Sym),nRoots)
  integer(kind=iwp) :: i, j
  real(kind=wp) :: scal
  real(kind=wp), external :: ddot_

  ! Project out the internal rotation contribution

  do i=1,nRoots
    do j=1,nRoots
      scal = DDot_(ncsf(State_Sym),CI(:,i),1,CIDER(:,j),1)
      CIDER(:,j) = CIDER(:,j)-Scal*CI(:,i)
    end do
  end do

end subroutine ISR_projection

!-----------------------------------------------------------------------

subroutine ISR_TimesE2(MODE,CI,CIDER)

  use input_mclr, only: ERASSCF, ncsf, nRoots, State_Sym, Weight
  !use DWSol, only: DWSCF, DWSol_Der

  integer(kind=iwp), intent(in) :: MODE
  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots)
  real(kind=wp), intent(inout) :: CIDER(ncsf(State_Sym),nRoots)
  integer(kind=iwp) :: i, j
  real(kind=wp) :: fact, scal
  !real(kind=wp), allocatable :: DERHII(:), DEROMG(:)
  real(kind=wp), external :: ddot_

# include "macros.fh"
  unused_var(mode)

  ! Compute the C-S and S-S blocks of the Ap operation

  !! if do_RF and unequal_SA are true, this subroutine is called twice, so the S-S block has to be halved
  fact = One
  if (do_RF .and. unequal_SA) fact = Half

  if (.not. InvSCF) then
    do i=1,nRoots
      do j=1,i-1
        scal = Zero
        if (IntRotOff) then
          !! Note that CIDER has been multiplied by Weight or W_SOLV
          scal = DDot_(ncsf(State_Sym),CI(:,i),1,CIDER(:,j),1)-DDot_(ncsf(State_Sym),CI(:,j),1,CIDER(:,i),1)
        end if
        if (ScalWeight .and. (abs(Weight(i)-Weight(j)) > 1.0e-9_wp)) then
          ISR%Ap(i,j) = ISR%Ap(i,j)+scal+(ERASSCF(i)-ERASSCF(j))*ISR%p(i,j)*Two*fact*(Weight(i)-Weight(j))
        else
          ISR%Ap(i,j) = ISR%Ap(i,j)+scal+(ERASSCF(i)-ERASSCF(j))*ISR%p(i,j)*Two*fact
        end if
      end do
      !! The diagonal contribution is for dynamically weighted methods
      ISR%Ap(i,i) = ISR%Ap(i,i)+ISR%p(i,i)*fact
    end do
  end if

  !! Derivative of H for DW-MCSCF is evaluated with CI derivatives
  !! DW solvation is evaluated in DWder_MCLR
  !if ((mode == 1) .and. DWSCF%do_DW) then
  !  call mma_allocate(DEROMG,nRoots,Label='DEROMG')
  !  DEROMG = Zero
  !  if (IntRotOff) then
  !    do i=1,nRoots
  !      ! CIDER has been scaled with the weight in cisigma_sa
  !      if (weight(i) >= 1.0e-8_wp) then
  !        DEROMG(i) = DDot_(ncsf(State_Sym),CI(:,i),1,CIDER(:,i),1)/weight(i)
  !      else
  !        DEROMG(i) = Zero !! under consideration
  !      end if
  !    end do
  !  end if
  !
  !  call mma_allocate(DERHII,nRoots,Label='DERHII')
  !  DERHII(:) = Zero
  !  call DWSol_Der(mode,DEROMG,DERHII,ERASSCF,weight)
  !
  !  do i=1,nRoots
  !    !! not sure why 1/4
  !    !! One reason is CIDER has been scaled by two, and the other is ?
  !    ISR%Ap(i,i) = ISR%Ap(i,i)+DERHII(i)*Quart
  !  end do
  !  call mma_deallocate(DERHII)
  !  call mma_deallocate(DEROMG)
  !end if

end subroutine ISR_TimesE2

!-----------------------------------------------------------------------

subroutine DMInvISR(ISRotIn,ISRotOut)

  use input_mclr, only: ERASSCF, nRoots, Weight

  real(kind=wp), intent(in) :: ISRotIn(nRoots,nRoots)
  real(kind=wp), intent(inout) :: ISRotOut(nRoots,nRoots)
  integer(kind=iwp) :: i, j
  logical(kind=iwp) :: Edeg, Wdeg

  ! diagonal preconditioning for the internal state rotation Hessian

  do i=1,nRoots
    do j=1,i-1
      Edeg = .false.
      if (abs(ERASSCF(i)-ERASSCF(j)) < 1.0e-9_wp) Edeg = .true.
      Wdeg = .false.
      if (ScalWeight .and. (abs(Weight(i)-Weight(j)) < 1.0e-9_wp)) Wdeg = .true.
      ! How to avoid the zero divsion for degenerated states?
      if (Edeg .and. Wdeg) then
        ISRotOut(i,j) = Half*ISRotIn(i,j)
      else if (Wdeg .or. (.not. ScalWeight)) then
        ISRotOut(i,j) = Half*ISRotIn(i,j)/(ERASSCF(i)-ERASSCF(j))
      else if (Edeg .and. ScalWeight) then
        ISRotOut(i,j) = Half*ISRotIn(i,j)/(Weight(i)-Weight(j))
      else
        ISRotOut(i,j) = Half*ISRotIn(i,j)/((ERASSCF(i)-ERASSCF(j))*(Weight(i)-Weight(j)))
      end if
    end do
    ISRotOut(i,i) = ISRotIn(i,i)
  end do

end subroutine DMInvISR

end module ISRotation
