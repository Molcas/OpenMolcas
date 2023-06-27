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
!***********************************************************************

subroutine RS_I_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc,Thr_RS)
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
!***********************************************************************

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
real*8 H(nInter,nInter), g(nInter), dq(nInter)
real*8, allocatable :: Val(:), Tmp(:,:), Vec(:,:), Mat(:)
character*6 UpMeth
character*1 Step_Trunc
logical Found

Lu = 6
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt(' In RS_I_RFO: H','(10f10.6)',H,nInter,nInter)
call RecPrt(' In RS_I_RFO: g','(10f10.6)',g,nInter,1)
call RecPrt(' In RS_I_RFO:dq','(10f10.6)',dq,nInter,1)
#endif

NumVal = min(2,nInter)
nVStep = 2
Found = .false.
Thr = 1.0D-6
call mma_allocate(Vec,nInter,NumVal,Label='Vec')
Vec(:,:) = Zero
call mma_allocate(Val,NumVal,Label='Val')
Val(:) = Zero
call mma_allocate(Mat,nInter*(nInter+1)/2,Label='Mat')
do i=1,nInter
  do j=1,i
    ij = i*(i-1)/2+j
    Mat(ij) = H(i,j)
  end do
end do

! Find the negative eigenvalue(s)
! Stop when the highest eigenvalue found is larger than Thr
do while (.not. Found)
  call Davidson(Mat,nInter,NumVal,Val,Vec,iStatus)
  if (iStatus > 0) then
    call SysWarnMsg('RS_I_RFO','Davidson procedure did not converge','')
  end if
  nNeg = 0
  do i=1,NumVal
    if (Val(i) < Zero) nNeg = nNeg+1
  end do
  if (((Val(NumVal) > Thr) .and. (nNeg > 0)) .or. (NumVal >= nInter)) then
    Found = .true.
  else
    ! Increase the number of eigenpairs to compute
    call mma_allocate(Tmp,nInter,NumVal,Label='Tmp')
    call dcopy_(NumVal*nInter,Vec,1,Tmp,1)
    call mma_deallocate(Vec)
    call mma_deallocate(Val)

    NumVal = min(NumVal+nVStep,nInter)

    call mma_allocate(Vec,nInter,NumVal,Label='Vec')
    call mma_allocate(Val,NumVal,Label='Val')
    Vec(:,:) = Zero
    Vec(:,1:NumVal-nVStep) = Tmp(:,:)
    Val(:) = Zero
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
if (nNeg == 0) then
  write(Lu,*) 'Warning RS-I-RFO: Neq == 0'
  call RecPrt(' In RS_I_RFO: Eigenvalues',' ',Val,1,NumVal)
  !call Abend()
end if
#ifdef _DEBUGPRINT_
call RecPrt(' In RS_I_RFO: Eigenvalues',' ',Val,1,NumVal)
call RecPrt(' In RS_I_RFO: Eigenvectors',' ',Vec,nInter,NumVal)
write(Lu,*) ' nNeg=',nNeg
#endif

! Transform the gradient and Hessian to generate the
! corresponding entities for the image function. This
! corresponds to an elementary Householder orthogonal
! transformation.

if (nNeg > 0) then
  call mma_allocate(Tmp,nInter,1,Label='Tmp')
  call dcopy_(nInter,g,1,Tmp(1,1),1)

  do iNeg=1,nNeg
    gi = DDot_(nInter,g,1,Vec(1,iNeg),1)
    call DaXpY_(nInter,-Two*gi,Vec(1,iNeg),1,g,1)
    Fact = Two*Val(iNeg)
    do j=1,nInter
      do i=1,nInter
        H(i,j) = H(i,j)-Fact*Vec(i,iNeg)*Vec(j,iNeg)
      end do
    end do
  end do
end if

call mma_deallocate(Vec)
call mma_deallocate(Val)

call RS_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,Step_Trunc,Thr_RS)

! Restore the original gradient

if (nNeg > 0) then
  call dcopy_(nInter,Tmp(1,1),1,g,1)
  call mma_deallocate(Tmp)
end if

UpMeth = 'RSIRFO'

#ifdef _DEBUGPRINT_
call RecPrt(' In RS_I_RFO: g','(10f10.6)',g,nInter,1)
call RecPrt(' In RS_I_RFO:dq','(10f10.6)',dq,nInter,1)
#endif

return

end subroutine RS_I_RFO
