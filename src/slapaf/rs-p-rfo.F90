!***********************************************************************
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
! Copyright (C) 1994,1997, Roland Lindh                                *
!               2014, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine RS_P_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc)
!***********************************************************************
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December '94                                             *
!                                                                      *
!                                                                      *
!             Solve      | H    g | | d |     | 1 0 | | d |            *
!                        |  T     | |   | = e |     | |   |            *
!                        | g    0 | | 1 |     | 0 1 | | 1 |            *
!                                                                      *
!             this corresponds to                                      *
!                                                                      *
!             H d + g = e d                                            *
!                                                                      *
!             and                                                      *
!                                                                      *
!               T                                                      *
!             g  d = e                                                 *
!                                                                      *
!             Modified from single negative eigenvalue to an arbitrary *
!             number, June '97, R. Lindh                               *
!                                                                      *
!             Removed full diagonalizations, April '14, I. Fdez. Galvan*
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nInter
real(kind=wp) :: H(nInter,nInter), g(nInter), dq(nInter), dqHdq, StepMax
character(len=6) :: UpMeth
character :: Step_Trunc
#include "print.fh"
integer(kind=iwp) :: i, ij, iPrint, iRout, iStatus, Iter, IterMx, j, jk, k, Lu, mInter, nNeg, NumVal, nVStep
real(kind=wp) :: A_RFO, A_RFO_long, A_RFO_short, dqdq, dqdq_long, dqdq_short, EigVal_r, EigVal_t, gv, Lambda, Thr
logical(kind=iwp) :: Found, Iterate
real(kind=wp), allocatable :: GradN(:), GradP(:), Mat(:), MatN(:), MatP(:), StepN(:), StepP(:), Tmp(:,:), TmpN(:), TmpP(:), &
                              Val(:), ValN(:), ValP(:), Vec(:,:), VecN(:), VecP(:)
real(kind=wp), external :: DDot_

iRout = 215
Lu = u6
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In RS_P_RFO: H','(10f10.6)',H,nInter,nInter)
  call RecPrt(' In RS_P_RFO: g','(10f10.6)',g,nInter,1)
  call RecPrt(' In RS_P_RFO:dq','(10f10.6)',dq,nInter,1)
end if

UpMeth = 'RSPRFO'

NumVal = min(2,nInter)
nVStep = 2
Found = .false.
Thr = 1.0e-6_wp
call mma_allocate(Vec,nInter,NumVal,Label='Vec')
call mma_allocate(Val,NumVal,Label='Val')
call mma_allocate(Mat,nInter*(nInter+1)/2,Label='Mat')
Vec(:,:) = Zero
do i=1,nInter
  do j=1,i
    ij = i*(i-1)/2+j
    Mat(ij) = H(i,j)
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the negative eigenvalue(s)
! Stop when the highest eigenvalue found is larger than Thr
do while (.not. Found)
  call Davidson(Mat,nInter,NumVal,Val,Vec,iStatus)
  if (iStatus > 0) then
    call SysWarnMsg('RS_P_RFO','Davidson procedure did not converge','')
  end if
  if ((Val(NumVal) > Thr) .or. (NumVal >= nInter)) then
    Found = .true.
  else
    ! Increase the number of eigenpairs to compute
    call mma_allocate(Tmp,nInter,NumVal,Label='Tmp')
    Tmp(:,:) = Vec(:,:)
    call mma_deallocate(Vec)
    call mma_deallocate(Val)

    NumVal = min(NumVal+nVStep,nInter)

    call mma_allocate(Vec,nInter,NumVal,Label='Vec')
    Vec(:,:) = Zero
    call mma_allocate(Val,NumVal,Label='Val')
    Val(:) = Zero

    Vec(:,1:NumVal-nVStep) = Tmp(:,:)

    call mma_deallocate(Tmp)
  end if
end do
call mma_deallocate(Mat)

nNeg = 0
i = NumVal
do while ((i >= 0) .and. (nNeg == 0))
  if (Val(i) < Zero) nNeg = i
  i = i-1
end do
if (iPrint >= 99) then
  call RecPrt(' In RS_P_RFO: Eigenvalues',' ',Val,1,NumVal)
  call RecPrt(' In RS_P_RFO: Eigenvectors',' ',Vec,nInter,NumVal)
  write(Lu,*) ' nNeg=',nNeg
end if

if (iPrint >= 6) then
  write(Lu,*)
  write(Lu,*) 'RS-P-RF Optimization'
  write(Lu,*) ' Iter   alpha        dqdq  StepMax   EigVal_r  EigVal_t'
end if

!write(Lu,*) 'Trust radius=',StepMax
A_RFO = One   ! Initial seed of alpha
IterMx = 25
Iter = 0
Iterate = .false.
Thr = 1.0e-7_wp
if (nNeg > 0) then
  mInter = nNeg+1
  call mma_allocate(StepN,nInter,Label='StepN')
  call mma_allocate(GradN,nInter,Label='GradN')
  call mma_allocate(VecN,mInter,Label='VecN')
  call mma_allocate(ValN,1,Label='ValN')
  call mma_allocate(MatN,mInter*(mInter+1)/2,Label='MatN')
  call mma_allocate(TmpN,mInter,Label='TmpN')
  TmpN(:) = Zero
end if
mInter = nInter+1
call mma_allocate(StepP,nInter,Label='StepP')
call mma_allocate(GradP,nInter,Label='GradP')
call mma_allocate(VecP,mInter,Label='VecP')
call mma_allocate(ValP,1,Label='ValP')
call mma_allocate(MatP,mInter*(mInter+1)/2,Label='MatP')
call mma_allocate(TmpP,mInter,Label='TmpP')
TmpP(:) = Zero
998 continue
Iter = Iter+1
!write(Lu,*) 'Iter=',Iter
!write(Lu,*) 'A_RFO=',A_RFO
call FZero(dq,nInter)
if (nNeg > 0) then
  !write(Lu,*)
  !write(Lu,*) 'Process negative eigenvalues.'
  !write(Lu,*)
  mInter = nNeg+1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Build the augmented matrix of the "negative" subspace
  !  diagonal: negative eigenvalues (divided by alpha)
  !  last row/column: components of the gradient along the
  !                   negative eigenvectors (divided by sqrt(alpha))
  ! Project the gradient in the negative subspace (but expressed
  ! in the full space)
  ! (Note that, since we are interested in the largest eigenvalue,
  !  the whole matrix is multiplied by -1, so we actually find the
  !  smallest eigenvalue and change its sign)
  GradN(:) = Zero
  MatN(:) = Zero
  j = mInter*(mInter-1)/2
  do i=1,nNeg
    MatN(i*(i+1)/2) = -Val(i)/A_RFO
    gv = DDot_(nInter,g,1,Vec(:,i),1)
    MatN(j+i) = gv/sqrt(A_RFO)
    call DaXpY_(nInter,gv,Vec(:,i),1,GradN,1)
  end do

  ! Solve the partial RFO system for the negative subspace
  VecN(:) = TmpN(:)
  call Davidson(MatN,mInter,1,ValN,VecN,iStatus)
  TmpN(:) = VecN(:)
  if (iStatus > 0) then
    call SysWarnMsg('RS_P_RFO','Davidson procedure did not converge','')
  end if
  ValN(1) = -ValN(1)

  ! Scale the eigenvector (combines eqs. (5) and (23))
  ! Convert to full space and add to complete step
  call DScal_(nNeg,One/(sqrt(A_RFO)*VecN(1+nNeg)),VecN,1)
  call dGeMV_('N',nInter,nNeg,One,Vec,nInter,VecN,1,Zero,StepN,1)
  call DaXpY_(nInter,One,StepN,1,dq,1)
  !dqdq_max = sqrt(DDot_(nInter,StepN,1,StepN,1))
  !write(Lu,*) 'dqdq_max=',dqdq_max
  ! Sign
  EigVal_r = -DDot_(nInter,StepN,1,GradN,1)
  if (iPrint >= 99) then
    call RecPrt('dq_r',' ',StepN,1,nInter)
    call RecPrt(' g_r',' ',GradN,1,nInter)
    write(Lu,*) 'Lambda=',EigVal_r
  end if
  if (EigVal_r < -Thr) then
    write(Lu,*)
    write(Lu,*) 'W A R N I N G !'
    write(Lu,*) 'EigVal_r < Zero',EigVal_r
    write(Lu,*)
  end if
else
  EigVal_r = Zero
  !dqdq_max = Zero
end if

!write(Lu,*)
!write(Lu,*) 'Process positive eigenvalues.'
!write(Lu,*)
mInter = nInter+1
!                                                                      *
!***********************************************************************
!                                                                      *
! Build the augmented matrix of the "positive" subspace
! Instead of reducing the dimensions, the negative eigenvectors
! are simply projected out from the gradient and the eigenvalues
! are shifted to a large positive arbitrary value (10), to avoid
! interferences
GradP(:) = g(:)
do i=1,nNeg
  gv = DDot_(nInter,GradP(:),1,Vec(:,i),1)
  call DaXpY_(nInter,-gv,Vec(:,i),1,GradP(:),1)
end do
do j=1,nInter
  call dcopy_(j,H(1,j),1,MatP(1+j*(j-1)/2),1)
  do i=1,nNeg
    do k=1,j
      jk = j*(j-1)/2+k
      MatP(jk) = MatP(jk)-(Val(i)-Ten)*Vec(j,i)*Vec(k,i)
    end do
  end do
  call DScal_(j,One/A_RFO,MatP(1+j*(j-1)/2),1)
end do
call FZero(MatP(1+mInter*(mInter-1)/2),mInter)
call DaXpY_(nInter,-One/sqrt(A_RFO),GradP(:),1,MatP(1+mInter*(mInter-1)/2),1)

! Solve the partial RFO system for the positive subspace
call dcopy_(mInter,TmpP(:),1,VecP(:),1)
call Davidson(MatP,mInter,1,ValP,VecP,iStatus)
if (iStatus > 0) then
  call SysWarnMsg('RS_P_RFO','Davidson procedure did not converge','')
end if
TmpP(1:mInter) = VecP(1:mInter)
StepP(1:nInter) = VecP(1:nInter)

! Scale the eigenvector (combines eqs. (5) and (23))
! Add to complete step
call DScal_(nInter,One/(sqrt(A_RFO)*VecP(1+nInter)),StepP,1)
call DaXpY_(nInter,One,StepP,1,dq,1)
! dqdq_min = sqrt(DDot_(nInter,StepP,1,StepP,1))
!write(Lu,*) 'dqdq_min=',dqdq_min
EigVal_t = -DDot_(nInter,StepP,1,GradP,1) ! Sign
if (iPrint >= 99) then
  call RecPrt('dq_t',' ',StepP,1,nInter)
  call RecPrt(' g_t',' ',GradP,1,nInter)
  write(Lu,*) 'Lambda=',EigVal_t
end if
if (EigVal_t > Thr) then
  write(Lu,*)
  write(Lu,*) 'W A R N I N G !'
  write(Lu,*) 'EigVal_t > Zero',EigVal_t
  write(Lu,*)
end if

Lambda = EigVal_t+EigVal_r
dqdq = sqrt(DDot_(nInter,dq,1,dq,1))

if (iPrint >= 6) write(Lu,'(I5,5F10.5)') Iter,A_RFO,dqdq,StepMax,EigVal_r,EigVal_t
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize data for iterative scheme (only at first iteration)

if (.not. Iterate) then
  A_RFO_long = A_RFO
  dqdq_long = dqdq
  A_RFO_short = Zero
  dqdq_short = dqdq_long+One
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! RF with constraints. Start iteration scheme if computed step is too long.

if ((Iter == 1) .and. (dqdq > StepMax)) Iterate = .true.
!                                                                      *
!***********************************************************************
!                                                                      *
! Procedure if the step length is not equal to the trust radius

if (Iterate .and. (abs(StepMax-dqdq) > Thr)) then
  Step_Trunc = '*'
  !write(Lu,*) 'StepMax-dqdq=',StepMax-dqdq
  call Find_RFO_Root(A_RFO_long,dqdq_long,A_RFO_short,dqdq_short,A_RFO,dqdq,StepMax)
  if (Iter > IterMx) then
    write(Lu,*) ' Too many iterations in RF'
    Go To 997
  end if
  Go To 998
end if

997 continue

call mma_deallocate(Vec)
call mma_deallocate(Val)

if (nNeg > 0) then
  mInter = nNeg+1
  call mma_deallocate(StepN)
  call mma_deallocate(GradN)
  call mma_deallocate(VecN)
  call mma_deallocate(ValN)
  call mma_deallocate(MatN)
  call mma_deallocate(TmpN)
end if
mInter = nInter+1
call mma_deallocate(StepP)
call mma_deallocate(GradP)
call mma_deallocate(VecP)
call mma_deallocate(ValP)
call mma_deallocate(MatP)
call mma_deallocate(TmpP)

if (iPrint >= 6) then
  write(Lu,*)
  write(Lu,*)
  write(Lu,*) 'Rational Function Optimization: Lambda=',Lambda
  write(Lu,*)
end if
dqHdq = dqHdq+Lambda*Half

if (iPrint >= 99) then
  write(Lu,*) 'EigVal,dqHdq=',Lambda,dqHdq
  call RecPrt(' In RS_P_RFO: g','(10f10.6)',g,nInter,1)
  call RecPrt(' In RS_P_RFO:dq','(10f10.6)',dq,nInter,1)
end if

return

end subroutine RS_P_RFO
