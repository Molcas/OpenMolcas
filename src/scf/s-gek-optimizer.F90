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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _FULL_SPACE_ ! Debugging
!#define _KRYLOV_
!#define _HYBRID_
!#define _HYBRID2_
#define _HYBRID3_
subroutine S_GEK_Optimizer(dq,mOV,dqdq,UpMeth,Step_Trunc)
!***********************************************************************
!                                                                      *
!     Object: subspace gradient-enhanced kriging optimization.         *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemistry -- BMC,                 *
!             University of Uppsala, SWEDEN                            *
!             May '22, November-December '22                           *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use InfSCF, only: Energy, HDiag, iter, iterso
use LnkLst, only: Init_LLs, LLGrad, LLx, LstPtr, SCF_V
use Kriging_mod, only: blaAI, blAI, blavAI, mblAI
use Kriging_procedures, only: Setup_Kriging
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Six, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mOV
real(kind=wp), intent(inout) :: dq(mOV)
real(kind=wp), intent(out) :: dqdq
character(len=6), intent(inout) :: UpMeth
character, intent(inout) :: Step_Trunc
integer(kind=iwp) :: i, iFirst, ii, ipg, ipq, Iter_Save, Iteration, Iteration_Micro, Iteration_Total, IterSO_Save, j, k, l, mDIIS, &
                     nDIIS, nExplicit
real(kind=wp) :: Beta_Disp, dqHdq, FAbs, Fact, gg, RMS, RMSMx, StepMax, Variance(1)
real(kind=wp), allocatable :: aux_a(:), aux_b(:), dq_diis(:), e_diis(:,:), g(:,:), g_diis(:,:), H_Diis(:,:), HDiag_Diis(:), &
                              q(:,:), q_diis(:,:), Val(:), Vec(:,:)
logical(kind=iwp) :: Converged, Terminate
character(len=6) :: UpMeth_
character :: Step_Trunc_
integer(kind=iwp), parameter :: Max_Iter = 50, nWindow = 8
#ifdef _KRYLOV_
integer(kind=iwp), parameter :: nKrylov = 20
#endif
real(kind=wp), parameter :: Beta_Disp_Min = 5.0e-3_wp, Beta_Disp_Seed = 0.05_wp, StepMax_Seed = 0.1_wp, Thr_RS = 1.0e-7_wp, &
                            ThrGrd = 1.0e-7_wp
real(kind=wp), external :: DDot_

Beta_Disp = Beta_Disp_Seed
#ifdef _DEBUGPRINT_
write(u6,*) 'Enter S-GEK Optimizer'
#endif
if (.not. Init_LLs) then
  write(u6,*) 'Link list not initiated'
  call Abend()
end if

call mma_allocate(q,mOV,iterso,Label='q')
call mma_allocate(g,mOV,iterso,Label='g')

!Pick up coordinates and gradients in full space
iFirst = iter-min(iterso,nWindow)+1
j = 0
do i=iFirst,iter
  j = i-iFirst+1
  !write(u6,*) 'i,j,iter=',i,j,iter

  ! Coordinates
  ipq = LstPtr(i,LLx)
  q(:,j) = SCF_V(ipq)%A(:)

  ! Gradients
  ipg = LstPtr(i,LLGrad)
  g(:,j) = SCF_V(ipg)%A(:)

end do

nDIIS = iter-iFirst+1

if (nDIIS == 1) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'Exit S-GEK Optimizer'
# endif
  call mma_deallocate(g)
  call mma_deallocate(q)
  return
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'nWindow=',nWindow
write(u6,*) 'nDIIS=',nDIIS
write(u6,*) 'IterSO=',IterSO
call RecPrt('q',' ',q,mOV,nDIIS)
call RecPrt('g',' ',g,mOV,nDIIS)
call RecPrt('g(:,nDIIS)',' ',g(:,nDIIS),mOV,1)
#endif

#ifdef _FULL_SPACE_

! Set up the full space
nExplicit = mOV
call mma_allocate(e_diis,mOV,nExplicit,Label='e_diis')
e_diis(:,:) = Zero
do k=1,nExplicit
  e_diis(k,k) = One
end do

#elif defined (_KRYLOV_)

! Set up unit vectors corresponding to a Krylov subspace for Adx=g
nExplicit = min(nKrylov,mOV)
!nExplicit = mOV
call mma_allocate(e_diis,mOV,nExplicit,Label='e_diis')
call mma_allocate(Aux_a,mOV,Label='Aux_a')
call mma_allocate(Aux_b,mOV,Label='Aux_b')
j = 1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j+1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
do j=3,nExplicit
  call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
  e_diis(:,j) = Aux_b(:)/sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))
  Aux_a(:) = Aux_b(:)
end do
call mma_deallocate(Aux_b)
call mma_deallocate(Aux_a)

#elif defined (_HYBRID_)

! Set up unit vectors corresponding to the subspace which the BFGS update will span.
nExplicit = min(2*(nDIIS-1)+nKrylov,mOV)
call mma_allocate(e_diis,mOV,nExplicit,Label='e_diis')

call mma_allocate(Aux_a,mOV,Label='Aux_a')
call mma_allocate(Aux_b,mOV,Label='Aux_b')
IterSO_save = IterSO
Iter_save = Iter
Iter = iFirst
IterSO = 1

j = 0
do k=1,nDIIS-1
  j = j+1
  Aux_a(:) = g(:,k+1)-g(:,k)
  e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))

  j = j+1
  Aux_a(:) = q(:,k+1)-q(:,k)
  call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
  e_diis(:,j) = Aux_b(:)/sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))

  iter = iter+1
  iterSO = iterSO+1
end do
IterSO = IterSO_save
Iter = Iter_save

! Add some unit vectors correponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j+1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j+1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
do k=j+1,nExplicit
  call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
  e_diis(:,k) = Aux_b(:)/sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))
  Aux_a(:) = Aux_b(:)
end do
call mma_deallocate(Aux_b)
call mma_deallocate(Aux_a)

#elif defined (_HYBRID2_)

! Set up unit vectors corresponding to the subspace which the BFGS update will span.
nExplicit = min(2*nDIIS-1+nKrylov,mOV)
call mma_allocate(e_diis,mOV,nExplicit,Label='e_diis')

call mma_allocate(Aux_a,mOV,Label='Aux_a')
call mma_allocate(Aux_b,mOV,Label='Aux_b')
IterSO_save = IterSO
Iter_save = Iter
Iter = iFirst
IterSO = 1

j = 0
do k=1,nDIIS
  j = j+1
  Aux_a(:) = g(:,k)
  e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))

  if (j == 2*nDIIS-1) exit
  j = j+1
  Aux_a(:) = q(:,k)
  call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
  call RecPrt('Aux_b(:)',' ',Aux_b(:),1,mOV)
  e_diis(:,j) = Aux_b(:)/sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))

  Iter = min(Iter+1,Iter_Save)
  IterSO = min(IterSO+1,IterSO_Save)
end do
IterSO = IterSO_save
Iter = Iter_save

! Add some unit vectors correponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j+1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j+1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
Aux_a(:) = e_diis(:,j)
do k=j+1,nExplicit
  call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
  e_diis(:,k) = Aux_b(:)/sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))
  Aux_a(:) = Aux_b(:)
  Aux_a(:) = e_diis(:,k)
end do
call mma_deallocate(Aux_b)
call mma_deallocate(Aux_a)

#elif defined (_HYBRID3_)

!nExplicit = 2 * (nDIIS - 1) + mOV + 2
nExplicit = 2*(nDIIS-1)+2
call mma_allocate(e_diis,mOV,nExplicit,Label='e_diis')

call mma_allocate(Aux_a,mOV,Label='Aux_a')
call mma_allocate(Aux_b,mOV,Label='Aux_b')
IterSO_save = IterSO
Iter_save = Iter
Iter = iFirst
IterSO = 1

j = 0
do k=1,nDIIS-1
  j = j+1
  Aux_a(:) = g(:,k+1)-g(:,k)
  e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))

  j = j+1
  Aux_a(:) = q(:,k+1)-q(:,k)
  !call SOrUpV(Aux_a(:),mOV,Aux_b(:),'GRAD','BFGS')
  Aux_b(:) = Aux_a(:)
  e_diis(:,j) = Aux_b(:)/sqrt(DDot_(mOV,Aux_b(:),1,Aux_b(:),1))

  iter = iter+1
  iterSO = iterSO+1
end do
IterSO = IterSO_save
Iter = Iter_save
call mma_deallocate(Aux_b)

! Add some unit vectors correponding to the Krylov subspace algorithm, g, Ag, A^2g, ....
j = j+1
Aux_a(:) = g(:,nDIIS)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
j = j+1
Aux_a(:) = dq(:)
e_diis(:,j) = Aux_a(:)/sqrt(DDot_(mOV,Aux_a(:),1,Aux_a(:),1))
call mma_deallocate(Aux_a)

#endif

! Now orthogonalize all unit vectors
#ifdef _DEBUGPRINT_
if (allocated(e_diis)) call RecPrt('e_diis(unnorm)',' ',e_diis,mOV,nExplicit)
#endif
do l=1,2
  j = 1
  do i=2,nExplicit
    do k=1,j
      gg = DDot_(mOV,e_diis(:,i),1,e_diis(:,k),1)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'i,k,gg=',i,k,gg
#     endif
      e_diis(:,i) = e_diis(:,i)-gg*e_diis(:,k)
    end do
    gg = DDot_(mOV,e_diis(:,i),1,e_diis(:,i),1) ! renormalize
#   ifdef _DEBUGPRINT_
    write(u6,*) 'j,i,gg=',j,i,gg
#   endif
    if (gg > 1.0e-17_wp) then   ! Skip vector if linear dependent.
      j = j+1
      e_diis(:,j) = e_diis(:,i)/sqrt(gg)
    end if
  end do
end do
mDIIS = j
#ifdef _DEBUGPRINT_
write(u6,*) '      mOV:',mOV
#ifdef _HYBRID_
write(u6,*) 'nExplicit:',nExplicit,'=',2*(nDIIS-1),'+',nKrylov
#elif defined (_HYBRID2_)
write(u6,*) 'nExplicit:',nExplicit,'=',2*nDIIS-1,'+',nKrylov
#else
write(u6,*) 'nExplicit:',nExplicit
#endif
write(u6,*) 'IterSO   :',IterSO
write(u6,*) '    nDIIS:',nDIIS
write(u6,*) '    mDIIS:',mDIIS

write(u6,*) 'Check the ortonormality'
do i=1,mDIIS
  do j=1,i
    write(u6,*) i,j,DDot_(mOV,e_diis(:,i),1,e_diis(:,j),1)
  end do
  write(u6,*)
end do
if (allocated(e_diis)) call RecPrt('e_diis',' ',e_diis,mOV,mDIIS)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected displacement coordinates. Note that the displacements are relative to the last coordinate, q(:,nDIIS). !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mma_allocate(q_diis,mDIIS,nDIIS+Max_Iter,Label='q_diis')
q_diis(:,:) = Zero
do i=1,nDIIS
  do k=1,mDIIS
    gg = Zero
    do l=1,mOV
      gg = gg+(q(l,i)-q(l,nDIIS))*e_diis(l,k)
    end do
    q_diis(k,i) = gg
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computed the projected gradients !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mma_allocate(g_diis,mDIIS,nDIIS+Max_Iter,Label='g_diis')
g_diis(:,:) = Zero
do i=1,nDIIS
  do k=1,mDIIS
    gg = Zero
    do l=1,mOV
      gg = gg+g(l,i)*e_diis(l,k)
    end do
    g_diis(k,i) = gg
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('q_diis',' ',q_diis,mDIIS,nDIIS)
call RecPrt('g_diis',' ',g_diis,mDIIS,nDIIS)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Project the approximate Hessian to the subspace !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_allocate(H_diis,mDIIS,mDIIS,Label='H_diis')
call mma_allocate(HDiag_diis,mDIIS,Label='HDiag_diis')

do i=1,mDiis
  do j=1,mDiis
    gg = Zero
    do l=1,mOV
      gg = gg+e_diis(l,i)*HDiag(l)*e_diis(l,j)
    end do
    H_diis(i,j) = gg
  end do
  HDiag_Diis(i) = H_Diis(i,i)
end do
#ifdef _DEBUGPRINT_
call RecPrt('H_diis(HDiag)',' ',H_diis,mDIIS,mDIIS)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_allocate(dq_diis,mDiis,Label='dq_Diis')

!We need to set the bias

blavAI = Ten
call Setup_Kriging(nDiis,mDiis,q_diis,g_diis,Energy(iFirst),Hessian_HMF=H_diis)
!call Setup_Kriging(nDiis,mDiis,q_diis,g_diis,Energy(iFirst),HDiag=HDiag_diis)
if (.false.) write(u6,*) blAI,mblAI,blaAI,blavAI
call mma_deallocate(HDiag_diis)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here starts the code doing the actual optimization

UpMeth = 'RVO'
Terminate = .false.
Step_Trunc = 'N'
Converged = .false.

Iteration = nDiis-1
Iteration_Micro = 0
Iteration_Total = iter-1
if (nDIIS > 1) Beta_Disp = min(Beta_Disp_Seed,max(Beta_Disp_Min,abs(Energy(iter)-Energy(iter-1))))
#ifdef _DEBUGPRINT_
write(u6,*) 'Energy(iter)-Energy(iter-1)=',Energy(iter)-Energy(iter-1)
write(u6,*) 'nDIIS=',nDIIS
write(u6,*) 'Beta_Disp_Seed=',Beta_Disp_Seed
write(u6,*) 'Beta_Disp_Min=',Beta_Disp_Min
write(u6,*) 'Beta_Disp=',Beta_Disp
#endif

do while ((.not. Converged) .and. (nDIIS > 1)) ! Micro iterate on the surrogate model

  Iteration_Micro = Iteration_Micro+1
  Iteration_Total = Iteration_Total+1
  Iteration = Iteration+1
  if (Iteration_Micro == Max_Iter) then
    write(u6,*)
    write(u6,*) 'S-GEK-Optimizer: Iteration_Micro==Max_Iter'
    write(u6,*) 'Abend!'
    write(u6,*)
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) '================================'
  write(u6,*) 'Micro Iteration=',Iteration_Micro
  write(u6,*) '================================'
  write(u6,*)
  write(u6,*) 'Step_Trunc:',Step_Trunc
  write(u6,*) '-----> Start RVO step'
# endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Fact = One
  StepMax = StepMax_Seed*real(Iteration_Micro,kind=wp)
  ! Loop to enforce restricted variance. Note, if the step restriction kicks in no problem since we will still microiterate.
  ! Normally a full step will be allowed -- no step restriction -- and the loop will be exited after the first iteration.
  do

    ! Compute the surrogate Hessian
    call Hessian_Kriging_Layer(q_diis(:,Iteration),H_diis,mDiis)
    !call Hessian_Kriging(q_diis(:,Iteration),H_diis,mDiis)

    call mma_allocate(Val,nTri_Elem(mDIIS),Label='Val')
    call mma_allocate(Vec,mDIIS,mDIIS,Label='Vec')

    call unitmat(Vec,mDIIS)
    do i=1,mDIIS
      do j=1,i
        Val(iTri(i,j)) = H_diis(i,j)
      end do
    end do

    call NIDiag_new(Val,Vec,mDIIS,mDIIS)
    call Jacord(Val,Vec,mDIIS,mDIIS)

    ! If negative eigenvalues then correct and signal that the micro iterartions should be terminanted.
    do i=1,mDIIS
      ii = nTri_Elem(i)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Eigenvalue:',Val(ii)
#     endif
      if (Val(ii) < Zero) then
        Terminate = .true.
        do j=1,mDIIS
          do k=1,mDIIS
            H_Diis(j,k) = H_Diis(j,k)+Two*abs(Val(ii))*Vec(j,i)*Vec(k,i)
          end do
        end do
      end if
    end do

    call mma_deallocate(Vec)
    call mma_deallocate(Val)

#   ifdef _DEBUGPRINT_
    call RecPrt('q_diis(:,Iteration)',' ',q_diis(:,Iteration),mDIIS,1)
    call RecPrt('H_diis(updated)',' ',H_diis,mDIIS,mDIIS)
    write(u6,*) 'Step_Trunc:',Step_Trunc
#   endif

    Step_Trunc_ = Step_Trunc
    dqHdq = Zero
    call RS_RFO(H_diis,g_Diis(:,Iteration),mDiis,dq_diis,UpMeth_,dqHdq,StepMax,Step_Trunc_,Thr_RS)
    dq_diis(:) = -dq_diis(:)
    q_diis(:,Iteration+1) = q_diis(:,Iteration)+dq_diis(:)
    dqdq = sqrt(DDot_(size(dq_diis),dq_diis(:),1,dq_diis(:),1))

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'Subiteration: Step_Trunc, StepMax,dqdq:',Step_Trunc,StepMax,dqdq
    write(u6,*) 'Subiteration: Step_Trunc_        :',Step_Trunc_
    call RecPrt('dq_diis',' ',dq_diis,mDIIS,1)
    call RecPrt('q_diis(:,Iteration+1)',' ',q_diis(:,Iteration+1),mDIIS,1)
#   endif
    if (Step_Trunc == 'N') Step_Trunc = ' '   ! set to blank if not touched
    if (Step_Trunc//Step_Trunc_ == ' *') Step_Trunc = '.' ! Mark that we have had a step Reduction

    call Dispersion_Kriging_Layer(q_diis(:,Iteration+1),Variance,mDIIS)
    !call Dispersion_Kriging(q_diis(:,Iteration+1),Variance,mDIIS)

    ! Note that we might have converged because the step restriction kicked in. However, we fill implicitly
    ! fix that during the second micro iteration.

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'Step_Trunc:',Step_Trunc
    write(u6,*) 'Beta_Disp =',Beta_Disp
    write(u6,*) 'Variance  =',Variance(1)
    write(u6,*) 'Fact      =',Fact
    write(u6,*) 'StepMax   =',StepMax
#   endif
    if (One-Variance(1)/Beta_Disp > 1.0e-3_wp) exit
    if ((Fact < 1.0e-5_wp) .or. (Variance(1) < Beta_Disp)) exit
    Fact = Half*Fact
    StepMax = Half*StepMax
    Step_Trunc = '*' ! This will only happen if the variance restriction kicks in

  end do  ! Restricted variance step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Step_Trunc:',Step_Trunc
  write(u6,*) '-----> Exit RVO step'
  write(u6,*)
# endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Compute the energy and the gradient of the surrogate model

  call Energy_Kriging_Layer(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
  !call Energy_Kriging(q_diis(:,Iteration+1),Energy(Iteration_Total+1),mDIIS)
  call Gradient_Kriging_Layer(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)
  !call Gradient_Kriging(q_diis(:,Iteration+1),g_diis(:,Iteration+1),mDIIS)

  !dEner = Energy(Iteration_Total+1)-Energy(Iteration_Total)

# ifdef _DEBUGPRINT_
  write(u6,*) 'Energy(Iteration_Total+1):',Energy(Iteration_Total+1)
  call RecPrt('g_diis(:,Iteration+1)',' ',g_diis(:,Iteration+1),mDIIS,1)
# endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Check on convergence criteria.

  FAbs = sqrt(DDot_(mDIIS,g_diis(:,Iteration+1),1,g_diis(:,Iteration+1),1)/real(mDIIS,kind=wp))
  RMS = sqrt(DDot_(mDIIS,dq_diis(:),1,dq_diis(:),1)/real(mDIIS,kind=wp))
  RMSMx = Zero
  do i=1,mDIIS
    RMSMx = max(RMSMx,abs(dq_diis(i)))
  end do
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'FAbs=',FAbs
  write(u6,*) 'RMS=',RMS
  write(u6,*) 'RMSMx=',RMSMx
  write(u6,*)
# endif
  if (Step_Trunc == '.') Step_Trunc = ' '
  Converged = ((FAbs < ThrGrd) .and. (RMS < Four*ThrGrd) .and. (RMSMx < ThrGrd*Six) .and. &
               ((Step_Trunc == ' ') .or. (Step_Trunc == '#')))
# ifdef _DEBUGPRINT_
  write(u6,*) 'Step_Trunc:',Step_Trunc
  write(u6,*) 'Converged:',Converged
# endif
  if (Step_Trunc == '*') Converged = .true.
  if (Terminate) then
    Step_Trunc = '#'
    exit
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do  ! While not converged

#ifdef _DEBUGPRINT_
write(u6,*) 'Converged'
write(u6,*) 'Energy(Iteration_Total+1):',Energy(Iteration_Total+1)
#endif
write(UpMeth(5:6),'(I2)') Iteration_Micro

! Compute the displacement in the reduced space relative to the last structure of the full space
dq_diis(:) = q_diis(:,Iteration+1)-q_diis(:,nDIIS)
! Compute the displacement in the full space.
dq(:) = Zero
do i=1,mDIIS
  dq(:) = dq(:)+dq_diis(i)*e_diis(:,i)
end do
dqdq = sqrt(DDot_(size(dq),dq(:),1,dq(:),1))

#ifdef _DEBUGPRINT_
call RecPrt('dq_diis',' ',dq_diis(:),size(dq_diis),1)
write(u6,*) '||dq||=',sqrt(DDot_(size(dq),dq(:),1,dq(:),1))
call RecPrt('dq',' ',dq(:),size(dq),1)
call RecPrt('g_diis(:,Iteration+1)',' ',g_diis(:,Iteration+1),size(g_diis,1),1)
#endif

call Finish_Kriging()
call mma_deallocate(dq_diis)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_deallocate(h_diis)
call mma_deallocate(q_diis)
call mma_deallocate(g_diis)
call mma_deallocate(e_diis,safe='*')

call mma_deallocate(g)
call mma_deallocate(q)

#ifdef _DEBUGPRINT_
write(u6,*) 'Exit S-GEK Optimizer'
#endif

end subroutine S_GEK_Optimizer
