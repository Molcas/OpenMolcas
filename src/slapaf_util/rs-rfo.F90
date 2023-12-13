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
! Copyright (C) 1994,2004,2014,2017,2019,2020, Roland Lindh            *
!               2014,2018, Ignacio Fdez. Galvan                        *
!***********************************************************************

subroutine RS_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc,Thr_RS)
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
!     Remove references to work, Roland Lindh                          *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(in) :: H(nInter,nInter), g(nInter), StepMax, Thr_RS
real(kind=wp), intent(out) :: dq(nInter)
character(len=6), intent(out) :: UpMeth
real(kind=wp), intent(inout) :: dqHdq
character, intent(out) :: Step_Trunc
integer(kind=iwp) :: i, ij, iRoot, iStatus, Iter, IterMx, iVal, j, Lu, NumVal
real(kind=wp) :: A_RFO, A_RFO_long, A_RFO_short, Dist, dqdq, dqdq_long, dqdq_short, EigVal, Fact, VV, ZZ
logical(kind=iwp) :: Iterate, Restart
real(kind=wp), allocatable :: Matrix(:), Val(:), Vec(:,:), Tmp(:)
#define _CHECK_UPDATE_
#ifdef _CHECK_UPDATE_
real(kind=wp), parameter :: Thr_Check = 1.0e2_wp
#endif
real(kind=wp), external :: DDot_

UpMeth = 'RS-RFO'
Lu = u6
!#define _DEBUGPRINT_
!#define _DEBUG2_
#ifdef _DEBUGPRINT_
call RecPrt(' In RS_RFO: H',' ',H,nInter,nInter)
call RecPrt(' In RS_RFO: g',' ',g,nInter,1)
write(Lu,*) 'Trust radius=',StepMax

write(Lu,*)
write(Lu,*) 'RS-RF Optimization'
write(Lu,*) ' Iter   alpha          dqdq    StepMax     EigVal'
#endif

A_RFO = One   ! Initial seed of alpha
IterMx = 25
Iter = 0
Iterate = .false.
Restart = .false.
NumVal = min(6,nInter)+1
call mma_allocate(Vec,nInter+1,NumVal,Label='Vec')
call mma_allocate(Val,NumVal,Label='Val')
call mma_allocate(Matrix,nTri_Elem(nInter+1),Label='Matrix')
call mma_allocate(Tmp,nInter+1,Label='Tmp')

Vec(:,:) = Zero
Tmp(:) = Zero
do
  Iter = Iter+1
# ifdef _DEBUG2_
  write(Lu,*) 'Iter=',Iter
  write(Lu,*) 'A_RFO=',A_RFO
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !      Execute step 1 of page 266                                    *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Set up the augmented Hessian of Eq. (2)

  ! Assume that the S-matrix is trivially diagonal

  do i=1,nInter
    do j=1,i
      ij = nTri_Elem(i-1)+j
      Matrix(ij) = Half*(H(i,j)+H(j,i))/A_RFO
    end do
  end do
  j = nInter+1
  do i=1,nInter
    ij = nTri_Elem(j-1)+i
    Matrix(ij) = -g(i)/sqrt(A_RFO) ! note sign
  end do
  Matrix(nTri_Elem(j)) = Zero
# ifdef _DEBUG2_
  call TriPrt('R_Tri',' ',Matrix,nInter+1)
# endif

  ! Restore the vector from the previous iteration, if any
  Vec(:,1) = Tmp(:)
  call Davidson(Matrix,nInter+1,NumVal,Val,Vec,iStatus)
  if (iStatus > 0) call SysWarnMsg('RS_RFO','Davidson procedure did not converge','')

  ! Pick up the root which represents the shortest displacement.

# ifdef _DEBUGPRINT_
  call RecPrt('Val',' ',Val,1,NumVal)
  call RecPrt('Vec',' ',Vec,nInter+1,NumVal)
# endif
  iRoot = -1
  Dist = 1.0e99_wp
  do iVal=1,NumVal
    if (Vec(nInter+1,iVal) == Zero) cycle
    VV = DDot_(nInter,Vec(1,iVal),1,Vec(1,iVal),1)
    ZZ = VV/A_RFO+Vec(nInter+1,iVal)**2
    Fact = Vec(nInter+1,iVal)/sqrt(ZZ)
    dqdq = VV/(A_RFO*Fact**2*ZZ)
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'iVal,A_RFO=',iVal,A_RFO
    write(u6,*) 'ZZ=',ZZ
    write(u6,*) 'Fact=',Fact
    write(u6,*) 'dqdq=',dqdq
#   endif
    if (dqdq < Dist) then
      iRoot = iVal
      Dist = dqdq
    end if
  end do
  if (iRoot == -1) then
    write(u6,*)
    write(u6,*) 'RS-RFO: Illegal iroot value!'
    call Abend()
  end if
  Tmp(:) = Vec(:,iRoot)
  Vec(1:nInter,iRoot) = Vec(1:nInter,iRoot)/sqrt(A_RFO)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !      Execute step 2 on page 266                                    *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# ifdef _DEBUG2_
  write(Lu,*) ' RF eigenvalue=',Val
# endif
  ZZ = DDot_(nInter+1,Vec(1,iRoot),1,Vec(1,iRoot),1)
  Vec(:,iRoot) = Vec(:,iRoot)/sqrt(ZZ)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !      Execute step 3 of page 266                                    *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Copy v^k_{n,i}

  dq(:) = Vec(1:nInter,iRoot)

  ! Pick v^k_{1,i}

  Fact = Vec(nInter+1,iRoot)
# ifdef _DEBUG2_
  write(Lu,*) 'v^k_{1,i}=',Fact
# endif

  ! Normalize according to Eq. (5)

  dq(:) = dq(:)/Fact
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'iRoot=',iRoot
  write(u6,*) 'ZZ=',ZZ
  write(u6,*) 'Fact=',Fact
# endif

  ! Compute lambda_i according to Eq. (8a)

  EigVal = -DDot_(nInter,dq,1,g,1) ! note sign

  ! Compute R^2 according to Eq. (8c)

  dqdq = sqrt(DDot_(nInter,dq,1,dq,1))
# ifdef _DEBUGPRINT_
  write(Lu,'(I5,5(ES12.5,1x))') Iter,A_RFO,dqdq,StepMax,EigVal
  !write(Lu,*) 'StepMax-dqdq=',StepMax-dqdq
  !write(Lu,*) 'Thr_RS=',Thr_RS
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Initialize data for iterative scheme (only at first iteration)

  if ((.not. Iterate) .or. Restart) then
    A_RFO_long = A_RFO
    dqdq_long = dqdq
    A_RFO_short = Zero
    dqdq_short = dqdq_long+One
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! RF with constraints. Start iteration scheme if computed step is too long.

  if (((Iter == 1) .or. Restart) .and. (dqdq > StepMax)) then
    Iterate = .true.
    Restart = .false.
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Procedure if the step length is not equal to the trust radius

  if ((.not. Iterate) .or. (abs(StepMax-dqdq) <= Thr_RS)) exit
  Step_Trunc = '*'
# ifdef _DEBUG2_
  write(Lu,*) 'StepMax-dqdq=',StepMax-dqdq
# endif

  ! Converge if small interval

  if ((dqdq < StepMax) .and. (abs(A_RFO_long-A_RFO_short) < Thr_RS)) exit
  call Find_RFO_Root(A_RFO_long,dqdq_long,A_RFO_short,dqdq_short,A_RFO,dqdq,StepMax)
  if (A_RFO == -One) then
    !write(Lu,*) 'reset Step_Trunc'
    A_RFO = One
    Step_Trunc = ' '
    Restart = .true.
    Iterate = .false.
  end if
  if (Iter > IterMx) then
    write(Lu,*) ' Too many iterations in RF'
    exit
  end if
end do

call mma_deallocate(Tmp)
dqHdq = dqHdq+EigVal*Half
#ifdef _DEBUGPRINT_
write(Lu,*)
write(Lu,*) 'Rational Function Optimization, Lambda=',EigVal
write(Lu,*)
write(Lu,*) 'EigVal,dqHdq=',EigVal,dqHdq
call RecPrt(' In RS_RFO: g',' ',g,nInter,1)
call RecPrt(' In RS_RFO:dq',' ',dq,nInter,1)
#endif
#define _CHECK_UPDATE_
#ifdef _CHECK_UPDATE_
do i=1,nInter
  if (abs(dq(i)) > Thr_Check) then
    write(u6,*) 'RS_RFO: ABS(dq(i)) > Thr_Check'
    write(u6,*) '        Probably an error.'
    call Abend()
  end if
end do
#endif

call mma_deallocate(Vec)
call mma_deallocate(Val)
call mma_deallocate(Matrix)
!write(u6,*) 'dqdq=',dqdq,dqdq**2
!write(u6,*) 'StepMax=',StepMax,StepMax**2
!write(Lu,*) 'StepMax-dqdq=',StepMax-dqdq
!write(Lu,*) dqdq < StepMax

return

end subroutine RS_RFO
