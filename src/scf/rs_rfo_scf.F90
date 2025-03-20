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
! Copyright (C) 1994,2004,2014,2017, Roland Lindh                      *
!               2014,2018,2024, Ignacio Fdez. Galvan                   *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RS_RFO_SCF(g,nInter,dq,UpMeth,dqdq,dqHdq,StepMax_Seed,Step_Trunc,ValMin)
!***********************************************************************
!                                                                      *
!     Object: Automatic restricted-step rational functional            *
!             optimization.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December '94                                             *
!                                                                      *
!     Modified to the restricted-step RFO method of Besalu and Bofill. *
!     Ref: E. Besalu and J. M. Bofill, TCA, 100, 265-274 (1998), by    *
!     R. Lindh, Gyeongju, Korea.                                       *
!     Removed full diagonalizations, Ignacio Fdez. Galvan, Uppsala     *
!     Remove references to work, Roland Lindh, Harvard, Cambridge      *
!     Modified for SCF, Roland Lindh, Harvard, Cambridge               *
!***********************************************************************

use InfSCF, only: IterSO, kOptim
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One, Two, Three, Pi
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nInter, ValMin
real(kind=wp), intent(in) :: g(nInter), StepMax_Seed
real(kind=wp), intent(out) :: dq(nInter), dqdq
character(len=6), intent(out) :: UpMeth
real(kind=wp), intent(inout) :: dqHdq
character, intent(out) :: Step_Trunc
integer(kind=iwp) :: I, iRoot, iStatus, Iter_i, NumVal
real(kind=wp) :: A_RFO, A_RFO_Long, A_RFO_Short, DqDq_Long, DqDq_Short, EigVal, Fact, GG, Step_Lasttime = Pi, StepMax = One, Test, &
                 ZZ
logical(kind=iwp) :: Iterate, Restart, Unreliable
real(kind=wp), allocatable :: Tmp(:), Val(:), Vec(:,:)
integer(kind=iwp), parameter :: IterMx = 20
real(kind=wp), parameter :: Step_Factor = Three, StepMax_Min = 1.0e-2_wp, Thr = 1.0e-4_wp
real(kind=wp), external :: DDot_

UpMeth = 'RS-RFO'
Step_Trunc = ' '
gg = sqrt(DDot_(nInter,g,1,g,1))
#ifdef _DEBUGPRINT_
write(u6,*) 'StepMax_Seed=',StepMax_Seed
write(u6,*) 'Sqrt(gg)=',gg
#endif

if (Step_Lasttime == Pi) Step_Lasttime = Step_Lasttime/(gg*Step_Factor)

StepMax = min(Pi,StepMax_Seed*gg,Step_Lasttime*Step_Factor*gg)

! Make sure that step restriction is not too tight.
if (StepMax < StepMax_Min) StepMax = StepMax_Min
#ifdef _DEBUGPRINT_
write(u6,*) 'StepMax=',StepMax
#endif

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) '***************************************************'
write(u6,*) '********* S T A R T  O F  R S - R F O SCF *********'
write(u6,*) '***************************************************'
call NrmClc(g,nInter,'RS-RFO','g(n)')
write(u6,*) 'Trust radius=',StepMax

write(u6,*)
write(u6,*) 'RS-RF Optimization'
write(u6,*) ' Iter    alpha    Sqrt(dqdq)  StepMax    EigVal'
#endif

A_RFO = One   ! Initial seed of alpha
Iter_i = 0
Iterate = .false.
Restart = .false.
Unreliable = .false.
NumVal = min(ValMin,nInter+1)
!NumVal = nInter+1
call mma_allocate(Vec,nInter+1,NumVal,Label='Vec')
call mma_allocate(Val,NumVal,Label='Val')
call mma_allocate(Tmp,nInter+1,Label='Tmp')

Vec(:,:) = Zero
Tmp(:) = Zero
do
  Iter_i = Iter_i+1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !        Execute step 1 of page 266                                  *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Restore the vector from the previous iteration, if any
  Vec(:,1) = Tmp(:)

  ! Call special Davidson routine which do not require the
  ! augmented Hessian to be explicitly expressed but rather will
  ! handle the gradient and Hessian part separated. The gradient
  ! will be explicit, while the Hessian part will use an approach
  ! which computes Hc, where c is a trial vector, from an initial
  ! Hessian based on a diagonal approximation and a BFGS update.

  call Davidson_SCF(g,nInter,NumVal,A_RFO,Val,Vec,iStatus)
  if (iStatus > 0) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'rs_rfo_SCF: Davidson procedure did not converge.'
#   endif
    call resetBFGS()
  end if
  !write(u6,*) 'Val(:)=',Val(:)
  !write(u6,*) 'Vec(:,1)=',Vec(:,1)
  !write(u6,*) 'Vec(nInter+1,1)=',Vec(nInter+1,1)

  ! Select a root with a negative value close to the current point

  iRoot = 1
  dqdq = 1.0e10_wp
  do i=1,NumVal
    if (Val(i) < Zero) then
      Tmp(:) = Vec(:,i)/sqrt(A_RFO)
      ZZ = DDot_(nInter+1,Tmp(:),1,Tmp(:),1)
      Tmp(:) = Tmp(:)/sqrt(ZZ)
      Tmp(:nInter) = Tmp(:nInter)/Tmp(nInter+1)
      Test = DDot_(nInter,Tmp,1,Tmp,1)
      !write(u6,*) 'Test=',Test
      if (Test < dqdq) then
        iRoot = i
        dqdq = Test
      end if
    end if
  end do
  !write(u6,*) 'iRoot,dqdq=',iRoot,dqdq
  if (iRoot /= 1) Vec(:,1) = Vec(:,iRoot)
  Tmp(:) = Vec(:,1)
  Vec(1:nInter,1) = Vec(1:nInter,1)/sqrt(A_RFO)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !        Execute step 2 on page 266                                  *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !write(u6,*) ' RF eigenvalue=',Val
  ZZ = DDot_(nInter+1,Vec(:,1),1,Vec(:,1),1)
  Vec(:,1) = Vec(:,1)/sqrt(ZZ)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !       Execute step 3 of page 266                                   *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Copy v^k_{n,i}

  dq(:) = Vec(1:nInter,1)

  ! Pick v^k_{1,i}

  Fact = Vec(nInter+1,1)
  !write(u6,*) 'v^k_{1,i}=',Fact

  ! Normalize according to Eq. (5)

  dq(:) = dq(:)/Fact

  ! Compute lambda_i according to Eq. (8a)

  EigVal = -DDot_(nInter,dq,1,g,1) ! note sign

  ! Compute R^2 according to Eq. (8c)

  dqdq = DDot_(nInter,dq,1,dq,1)

# ifdef _DEBUGPRINT_
  write(u6,'(I5,4ES11.3)') Iter_i,A_RFO,sqrt(dqdq),StepMax,EigVal
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Initialize data for iterative scheme (only at first iteration)

  if ((.not. Iterate) .or. Restart) then
    A_RFO_long = A_RFO
    dqdq_long = sqrt(dqdq)
    A_RFO_short = Zero
    dqdq_short = dqdq_long+One
  end if
  !write(u6,*) 'dqdq_long=',dqdq_long
  !write(u6,*) 'dqdq_short=',dqdq_short
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! RF with constraints. Start iteration scheme if computed step
  ! is too long.

  if (((Iter_i == 1) .or. Restart) .and. (dqdq > StepMax**2)) then
    Iterate = .true.
    Restart = .false.
  end if

  if (sqrt(dqdq) > Two*max(StepMax,Pi)) then
    !if ((sqrt(dqdq) > Pi) .or. (sqrt(dqdq) > StepMax) .and. (kOptim > 1)) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'rs_rfo_SCF: Total displacement is too large.'
#   endif
    call resetBFGS()
    if (Restart) cycle
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Procedure if the step length is not equal to the trust radius

  Step_Trunc = '*'
  !write(u6,*) 'StepMax-Sqrt(dqdq)=',StepMax-Sqrt(dqdq)

  ! Converge if small interval

  if (dqdq < StepMax**2) then
    ! Converge if we can't trust this
    if (Unreliable) exit
    ! Converge if small interval
    if (abs(A_RFO_long-A_RFO_short) < Thr) exit
  end if

  call Find_RFO_Root(A_RFO_long,dqdq_long,A_RFO_short,dqdq_short,A_RFO,sqrt(dqdq),StepMax)
  !write(u6,*) 'A_RFO_Short,A_RFO_Long,A_RFO=',A_RFO_Short,A_RFO_Long,A_RFO
  !write(u6,*) 'dqdq_short,dqdq_long,StepMax=',dqdq_short,dqdq_long,StepMax
  if (A_RFO == -One) then
    A_RFO = One
    Step_Trunc = ' '
    Restart = .true.
    Iterate = .false.
  end if

  !if (abs(StepMax-sqrt(dqdq)) <= Thr) write(6,*) 'Converged'
  if ((.not. Iterate) .or. (abs(StepMax-sqrt(dqdq)) <= Thr)) exit
  if (Iter_i > IterMx) exit

end do

! Safety net: truncate if the step was still too large
if (sqrt(dqdq)-StepMax > Thr) then
# ifdef _DEBUGPRINT_
  write(u6,*) ' Too many iterations in RS'
# endif
  dq(:) = dq(:)/sqrt(dqdq)
  dqdq = DDot_(nInter,dq,1,dq,1)
end if

call mma_deallocate(Tmp)
dqHdq = dqHdq+EigVal*Half
Step_Lasttime = sqrt(dqdq)/gg
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Rational Function Optimization, Lambda=',EigVal
write(u6,*)
write(u6,*) 'EigVal,dqHdq=',EigVal,dqHdq
call NrmClc(g,nInter,'RS-RFO','g(n)')
call NrmClc(dq,nInter,'RS-RFO','dX(n)')
write(u6,*) '***************************************************'
write(u6,*) '************* E N D  O F  R S - R F O SCF *********'
write(u6,*) '***************************************************'
write(u6,*)
#endif
call mma_deallocate(Vec)
call mma_deallocate(Val)

contains

subroutine resetBFGS()

  if (kOptim /= 1) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Reset update depth in BFGS, redo the RS-RFO'
#   endif
    kOptim = 1
    IterSO = 1
    A_RFO = One
    Restart = .true.
    Iterate = .false.
  else
    ! Something bad is going on
    Unreliable = .true.
  end if

end subroutine

end subroutine RS_RFO_SCF
