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
!               2025, Lila Zapp                                        *
!***********************************************************************

subroutine GEK_Optimizer(mDiis,nDiis,Max_Iter,q_diis,g_diis,dq_diis,Energy,H_diis,dqdq,Step_Trunc,UpMeth,SORange)

use Index_Functions, only: iTri, nTri_Elem
use Kriging_mod, only: blavAI
use Kriging_procedures, only: Setup_Kriging
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Six, Ten, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: mDiis, nDiis, Max_Iter
real(kind=wp), intent(inout) :: q_diis(mDiis,nDiis+Max_Iter), g_diis(mDiis,nDiis+Max_Iter), Energy(nDiis+Max_Iter), &
                                H_diis(mDiis,mDiis)
real(kind=wp), intent(out) :: dq_diis(mDiis), dqdq
character, intent(out) :: Step_Trunc
character(len=6), intent(out) :: UpMeth
logical(kind=iwp), intent(in) :: SORange
integer(kind=iwp) :: i, ii, Iteration, Iteration_Micro, Iteration_Total, j, k
real(kind=wp) :: Beta_Disp, dqHdq, FAbs, Fact, RMS, RMSMx, SOFact, StepMax, Variance(1)
logical(kind=iwp) :: Converged, Terminate
character(len=6) :: UpMeth_
character :: Step_Trunc_
real(kind=wp), allocatable :: Val(:), Vec(:,:)
integer(kind=iwp), parameter :: nWindow = 20
real(kind=wp), parameter :: Beta_Disp_Min = 5.0e-3_wp, Beta_Disp_Seed = 0.05_wp, StepMax_Seed = 0.1_wp, Thr_RS = 1.0e-7_wp, &
                            ThrGrd = 1.0e-7_wp
real(kind=wp), external :: DDot_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We need to set the bias

blavAI = Ten
call Setup_Kriging(nDiis,mDiis,q_diis,g_diis,Energy,Hessian_HMF=H_diis)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here starts the code doing the actual optimization

UpMeth = 'RVO'
Terminate = .false.
Step_Trunc = 'N'
Converged = .false.

if (SORange) then
  SOFact = One
else
  SOFact = 10000.0_wp
end if
Beta_Disp = Beta_Disp_Seed*SOFact
Iteration = nDiis-1
Iteration_Micro = 0
Iteration_Total = nDIIS-1
if (nDIIS > 1) Beta_Disp = min(Beta_Disp_Seed*SOFact,max(Beta_Disp_Min,abs(Energy(nDIIS)-Energy(nDIIS-1))))
#ifdef _DEBUGPRINT_
write(u6,*) 'Energy(nDIIS)-Energy(nDIIS-1)=',Energy(nDIIS)-Energy(nDIIS-1)
write(u6,*) 'nDIIS=',nDIIS
write(u6,*) 'Beta_Disp_Seed=',Beta_Disp_Seed*SOFact
write(u6,*) 'Beta_Disp_Min=',Beta_Disp_Min
write(u6,*) 'Beta_Disp=',Beta_Disp
#endif

do while (.not. Converged) ! Micro iterate on the surrogate model

  Iteration_Micro = Iteration_Micro+1
  Iteration_Total = Iteration_Total+1
  Iteration = Iteration+1

  !if (Iteration_Micro == Max_Iter) then
  !  write(u6,*)
  !  write(u6,*) 'S_GEK_Optimizer: Iteration_Micro==Max_Iter'
  !  write(u6,*) 'Abend!'
  !  write(u6,*)
  !  call Abend()
  !end if

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
  StepMax = StepMax_Seed*SOFact*real(Iteration_Micro,kind=wp)
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

    ! If negative eigenvalues then correct and signal that the micro iterations should be terminanted.
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

    ! Note that we might have converged because the step restriction kicked in. However, we will implicitly
    ! fix that during the second micro iteration.

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'Step_Trunc:',Step_Trunc
    write(u6,*) 'Beta_Disp =',Beta_Disp
    write(u6,*) 'Variance  =',Variance(1)
    write(u6,*) 'Fact      =',Fact
    write(u6,*) 'StepMax   =',StepMax
#   endif
    if ((Fact < 1.0e-5_wp) .or. (Variance(1) < Beta_Disp)) exit
    if (One-Variance(1)/Beta_Disp > 1.0e-3_wp) exit
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
  if ((.not. Converged) .and. (Iteration_Micro == Max_Iter)) Terminate = .true.
  if (Terminate) then
    Step_Trunc = '#'
    exit
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do  ! While not converged

#ifdef _DEBUGPRINT_
if (Converged) then
  write(u6,*) 'Converged'
else
  write(u6,*) 'Not converged!'
end if
write(u6,*) 'Energy(Iteration_Total+1):',Energy(Iteration_Total+1)
#endif
write(UpMeth(5:6),'(I2)') Iteration_Micro

! Compute the displacement in the reduced space relative to the last structure of the full space
dq_diis(:) = q_diis(:,Iteration+1)-q_diis(:,nDIIS)

#ifdef _DEBUGPRINT_
call RecPrt('dq_diis',' ',dq_diis(:),size(dq_diis),1)
call RecPrt('g_diis(:,Iteration+1)',' ',g_diis(:,Iteration+1),size(g_diis,1),1)
#endif

end subroutine GEK_Optimizer
