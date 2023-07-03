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

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
integer nInter
real*8 H(nInter,nInter), g(nInter), dq(nInter)
character UpMeth*6, Step_Trunc*1
real*8 StepMax
! Local variables
real*8, dimension(:), allocatable :: Tmp, Val, Matrix
real*8, dimension(:,:), allocatable :: Vec
logical Iterate, Restart

UpMeth = 'RS-RFO'
Lu = 6
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
call mma_allocate(Vec,(nInter+1),NumVal,Label='Vec')
call mma_allocate(Val,NumVal,Label='Val')
call mma_allocate(Matrix,(nInter+1)*(nInter+2)/2,Label='Matrix')
call mma_allocate(Tmp,nInter+1,Label='Tmp')

Vec(:,:) = 0.0d0
Tmp(:) = 0.0
998 continue
Iter = Iter+1
#ifdef _DEBUG2_
write(Lu,*) 'Iter=',Iter
write(Lu,*) 'A_RFO=',A_RFO
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!        Execute step 1 of page 266                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up the augmented Hessian of Eq. (2)

! Assume that the S-matrix is trivially diagonal

do i=1,nInter
  do j=1,i
    ij = i*(i-1)/2+j
    Matrix(ij) = Half*(H(i,j)+H(j,i))/A_RFO
  end do
end do
j = nInter+1
do i=1,nInter
  ij = j*(j-1)/2+i
  Matrix(ij) = -g(i)/sqrt(A_RFO) ! note sign
end do
jj = j*(j+1)/2
Matrix(jj) = Zero
#ifdef _DEBUG2_
call TriPrt('R_Tri',' ',Matrix,nInter+1)
#endif

! Restore the vector from the previous iteration, if any
call dcopy_(nInter+1,Tmp(1),1,Vec(1,1),1)
call Davidson(Matrix,nInter+1,NumVal,Val,Vec,iStatus)
if (iStatus > 0) then
  call SysWarnMsg('RS_RFO','Davidson procedure did not converge','')
end if

! Pick up the root which represents the shortest displacement.

#ifdef _DEBUGPRINT_
call RecPrt('Val',' ',Val,1,NumVal)
call RecPrt('Vec',' ',Vec,nInter+1,NumVal)
#endif
iRoot = -1
Dist = 1.0d99
do iVal=1,NumVal
  if (Vec(nInter+1,iVal) == 0.0d0) cycle
  VV = DDot_(nInter,Vec(1,iVal),1,Vec(1,iVal),1)
  ZZ = VV/A_RFO+Vec(nInter+1,iVal)**2
  Fact = Vec(nInter+1,iVal)/sqrt(ZZ)
  dqdq = VV/(A_RFO*Fact**2*ZZ)
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'iVal,A_RFO=',iVal,A_RFO
  write(6,*) 'ZZ=',ZZ
  write(6,*) 'Fact=',Fact
  write(6,*) 'dqdq=',dqdq
# endif
  if (dqdq < Dist) then
    iRoot = iVal
    Dist = dqdq
  end if
end do
if (iRoot == -1) then
  write(6,*)
  write(6,*) 'RS-RFO: Illegal iroot value!'
  call Abend()
end if
call dcopy_(nInter+1,Vec(1,iRoot),1,Tmp,1)
call DScal_(nInter,One/sqrt(A_RFO),Vec(1,iRoot),1)
!                                                                      *
!***********************************************************************
!                                                                      *
!        Execute step 2 on page 266                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUG2_
write(Lu,*) ' RF eigenvalue=',Val
#endif
ZZ = DDot_(nInter+1,Vec(1,iRoot),1,Vec(1,iRoot),1)
call DScal_(nInter+1,One/sqrt(ZZ),Vec(1,iRoot),1)
!                                                                      *
!***********************************************************************
!                                                                      *
!        Execute step 3 of page 266                                    *
!                                                                      *
!***********************************************************************
!                                                                      *
! Copy v^k_{n,i}

call dcopy_(nInter,Vec(1,iRoot),1,dq,1)

! Pick v^k_{1,i}

Fact = Vec(nInter+1,iRoot)
#ifdef _DEBUG2_
write(Lu,*) 'v^k_{1,i}=',Fact
#endif

! Normalize according to Eq. (5)

call DScal_(nInter,One/Fact,dq,1)
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'iRoot=',iRoot
write(6,*) 'ZZ=',ZZ
write(6,*) 'Fact=',Fact
#endif

! Compute lambda_i according to Eq. (8a)

EigVal = -DDot_(nInter,dq,1,g,1) ! note sign

! Compute R^2 according to Eq. (8c)

dqdq = sqrt(DDot_(nInter,dq,1,dq,1))
#ifdef _DEBUGPRINT_
write(Lu,'(I5,5(E12.5,1x))') Iter,A_RFO,dqdq,StepMax,EigVal
!write(Lu,*) 'StepMax-dqdq=',StepMax-dqdq
!write(Lu,*) 'Thr_RS=',Thr_RS
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize data for iterative scheme (only at first iteration)

if ((.not. Iterate) .or. Restart) then
  A_RFO_long = A_RFO
  dqdq_long = dqdq
  A_RFO_short = Zero
  dqdq_short = dqdq_long+One
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! RF with constraints. Start iteration scheme if computed step is too long.

if (((Iter == 1) .or. Restart) .and. (dqdq > StepMax)) then
  Iterate = .true.
  Restart = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Procedure if the step length is not equal to the trust radius

if (Iterate .and. (abs(StepMax-dqdq) > Thr_RS)) then
  Step_Trunc = '*'
# ifdef _DEBUG2_
  write(Lu,*) 'StepMax-dqdq=',StepMax-dqdq
# endif

  ! Converge if small interval

  if ((dqdq < StepMax) .and. (abs(A_RFO_long-A_RFO_short) < Thr_RS)) Go To 997
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
    Go To 997
  end if
  Go To 998
end if

997 continue
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
Thr_Check = 1.0d2
do i=1,nInter
  if (abs(dq(i)) > Thr_Check) then
    write(6,*) 'RS_RFO: ABS(dq(i)) > Thr_Check'
    write(6,*) '        Probably an error.'
    call Abend()
  end if
end do
#endif

call mma_deallocate(Vec)
call mma_deallocate(Val)
call mma_deallocate(Matrix)
!write(6,*) 'dqdq=',dqdq,dqdq**2
!write(6,*) 'StepMax=',StepMax,StepMax**2
!write(Lu,*) 'StepMax-dqdq=',StepMax-dqdq
!write(Lu,*) dqdq < StepMax

return

end subroutine RS_RFO
