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
!               2014,2018, Ignacio Fdez. Galvan                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RS_RFO_SCF(g,nInter,dq,UpMeth,dqdq,dqHdq,StepMax_Seed,Step_Trunc)
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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, One, Three, Pi
use InfSCF, only: kOptim, Iter_x => Iter, Iter_Start
use InfSO, only: IterSO

implicit none
integer nInter
real*8 g(nInter), dq(nInter)
character UpMeth*6
real*8 dqdq, dqHdq, StepMax_Seed
character Step_Trunc*1
! Local variables
integer :: Lu, IterMx, NumVal, iStatus, iRoot, I, Iter
real*8 :: GG, A_RFO, ZZ, Test, Fact, EigVal
real*8 :: A_RFO_Long, A_RFO_Short
real*8 :: DqDq_Long, DqDq_Short
real*8, external :: DDot_
real*8, allocatable :: Tmp(:), Val(:), Vec(:,:)
logical Iterate, Restart
real*8, save :: StepMax = One
real*8, save :: Step_Lasttime = Pi
real*8, parameter :: Thr = 1.0D-4
real*8, parameter :: StepMax_Min = 1.0D-2
real*8, parameter :: Step_Factor = Three

UpMeth = 'RS-RFO'
Step_Trunc = ' '
Lu = 6
gg = sqrt(DDot_(nInter,g,1,g,1))
#ifdef _DEBUGPRINT_
write(6,*) 'StepMax_Seed=',StepMax_Seed
write(6,*) 'Sqrt(gg)=',gg
#endif

if (Step_Lasttime == Pi) Step_Lasttime = Step_Lasttime/(gg*Step_Factor)

StepMax = min(Pi,StepMax_Seed*gg,Step_Lasttime*Step_Factor*gg)

! Make sure that step restriction is not too tight.
if (StepMax < StepMax_Min) StepMax = StepMax_Min
#ifdef _DEBUGPRINT_
write(6,*) 'StepMax=',StepMax
#endif

#ifdef _DEBUGPRINT_
write(Lu,*)
write(Lu,*) '***************************************************'
write(Lu,*) '********* S T A R T  O F  R S - R F O SCF *********'
write(Lu,*) '***************************************************'
call NrmClc(g,nInter,'RS-RFO','g(n)')
write(Lu,*) 'Trust radius=',StepMax

write(Lu,*)
write(Lu,*) 'RS-RF Optimization'
write(Lu,*) ' Iter    alpha    Sqrt(dqdq)  StepMax    EigVal'
#endif

A_RFO = One   ! Initial seed of alpha
IterMx = 50
Iter = 0
Iterate = .false.
Restart = .false.
NumVal = min(3,nInter+1)
!NumVal = min(nInter+1,nInter+1)
call mma_allocate(Vec,(nInter+1),NumVal,Label='Vec')
call mma_allocate(Val,NumVal,Label='Val')
call mma_allocate(Tmp,nInter+1,Label='Tmp')

Vec(:,:) = Zero
Tmp(:) = Zero
998 continue
Iter = Iter+1
!                                                                      *
!***********************************************************************
!                                                                      *
!        Execute step 1 of page 266                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
! Restore the vector from the previous iteration, if any
call dcopy_(nInter+1,Tmp,1,Vec(:,1),1)

! Call special Davidson routine which do not require the
! augmented Hessian to be explicitly expressed but rather will
! handle the gradient and Hessian part separated. The gradient
! will be explicit, while the Hessian part will use an approach
! which computes Hc, where c is a trial vector, from an initial
! Hessian based on a diagonal approximation and a BFGS update.

call Davidson_SCF(g,nInter,NumVal,A_RFO,Val,Vec,iStatus)
if (iStatus > 0) call SysWarnMsg('RS_RFO SCF','Davidson procedure did not converge','')
!write(6,*) 'Val(:)=',Val(:)
!write(6,*) 'Vec(:,1)=',Vec(:,1)
!write(6,*) 'Vec(nInter+1,1)=',Vec(nInter+1,1)

! Select a root with a negative value close to the current point

iRoot = 1
dqdq = 1.0d10
do i=1,NumVal
  if (Val(i) < Zero) then
    Tmp(:) = Vec(:,i)/sqrt(A_RFO)
    ZZ = DDot_(nInter+1,Tmp(:),1,Tmp(:),1)
    Tmp(:) = Tmp(:)/sqrt(ZZ)
    Tmp(:nInter) = Tmp(:nInter)/Tmp(nInter+1)
    Test = DDot_(nInter,Tmp,1,Tmp,1)
    !write(6,*) 'Test=',Test
    if (Test < dqdq) then
      iRoot = i
      dqdq = Test
    end if
  end if
end do
!write(6,*) 'iRoot,dqdq=',iRoot,dqdq
if (iRoot /= 1) Vec(:,1) = Vec(:,iRoot)
call dcopy_(nInter+1,Vec(:,1),1,Tmp,1)
call DScal_(nInter,One/sqrt(A_RFO),Vec(:,1),1)
!                                                                      *
!***********************************************************************
!                                                                      *
!        Execute step 2 on page 266                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
!write(Lu,*) ' RF eigenvalue=',Val
ZZ = DDot_(nInter+1,Vec(:,1),1,Vec(:,1),1)
call DScal_(nInter+1,One/sqrt(ZZ),Vec(:,1),1)
!                                                                      *
!***********************************************************************
!                                                                      *
!       Execute step 3 of page 266                                     *
!                                                                      *
!***********************************************************************
!                                                                      *
! Copy v^k_{n,i}

call dcopy_(nInter,Vec(:,1),1,dq,1)

! Pick v^k_{1,i}

Fact = Vec(nInter+1,1)
!write(Lu,*) 'v^k_{1,i}=',Fact

! Normalize according to Eq. (5)

call DScal_(nInter,One/Fact,dq,1)

! Compute lambda_i according to Eq. (8a)

EigVal = -DDot_(nInter,dq,1,g,1) ! note sign

! Compute R^2 according to Eq. (8c)

dqdq = DDot_(nInter,dq,1,dq,1)

if (sqrt(dqdq) > Pi) then
  !if ((sqrt(dqdq) > Pi) .or. (sqrt(dqdq) > StepMax) .and. (kOptim > 1)) then
  if (kOptim /= 1) then
    write(Lu,*) 'rs_rfo_SCF: Total displacement is too large.'
    write(Lu,*) 'DD=',sqrt(dqdq)
    write(Lu,*) 'Reset update depth in BFGS, redo the RS-RFO'
    Iter = Iter-1
    kOptim = 1
    Iter_Start = Iter_x
    IterSO = 1
    Go To 998
  end if
end if
#ifdef _DEBUGPRINT_
write(Lu,'(I5,4ES11.3)') Iter,A_RFO,sqrt(dqdq),StepMax,EigVal
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize data for iterative scheme (only at first iteration)

if ((.not. Iterate) .or. Restart) then
  A_RFO_long = A_RFO
  dqdq_long = sqrt(dqdq)
  A_RFO_short = Zero
  dqdq_short = dqdq_long+One
end if
!write(6,*) 'dqdq_long=',dqdq_long
!write(6,*) 'dqdq_short=',dqdq_short
!                                                                      *
!***********************************************************************
!                                                                      *
! RF with constraints. Start iteration scheme if computed step
! is too long.

if (((Iter == 1) .or. Restart) .and. (dqdq > StepMax**2)) then
  Iterate = .true.
  Restart = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Procedure if the step length is not equal to the trust radius

if (Iterate .and. (abs(StepMax-sqrt(dqdq)) > Thr)) then
  Step_Trunc = '*'
  !write(Lu,*) 'StepMax-Sqrt(dqdq)=',StepMax-Sqrt(dqdq)

  ! Converge if small interval

  if ((dqdq < StepMax**2) .and. (abs(A_RFO_long-A_RFO_short) < Thr)) Go To 997
  call Find_RFO_Root(A_RFO_long,dqdq_long,A_RFO_short,dqdq_short,A_RFO,sqrt(dqdq),StepMax)
  !write(6,*) 'A_RFO_Short=',A_RFO_Short
  !write(6,*) 'A_RFO_Long=',A_RFO_Long
  !write(6,*) 'dqdq_long=',dqdq_long
  !write(6,*) 'dqdq_short=',dqdq_short
  if (A_RFO == -One) then
    A_RFO = One
    Step_Trunc = ' '
    Restart = .true.
    Iterate = .false.
  end if
  if (Iter > IterMx) then
    write(Lu,*) ' Too many iterations in RF'
    Go To 997
  end if
  Go To 998
end if

997 continue
call mma_deallocate(Tmp)
dqHdq = dqHdq+EigVal*Half
Step_Lasttime = sqrt(dqdq)/gg
#ifdef _DEBUGPRINT_
write(Lu,*)
write(Lu,*) 'Rational Function Optimization, Lambda=',EigVal
write(Lu,*)
write(Lu,*) 'EigVal,dqHdq=',EigVal,dqHdq
call NrmClc(g,nInter,'RS-RFO','g(n)')
call NrmClc(dq,nInter,'RS-RFO','dX(n)')
write(Lu,*) '***************************************************'
write(Lu,*) '************* E N D  O F  R S - R F O SCF *********'
write(Lu,*) '***************************************************'
write(Lu,*)
#endif
call mma_deallocate(Vec)
call mma_deallocate(Val)

end subroutine RS_RFO_SCF
