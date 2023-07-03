!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ElRed(Bmtrx,nq,nx,Gmtrx,EVal,EVec,nK,uMtrx,Scrt,g12K,Thr)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
real*8 Bmtrx(nq,nx), Gmtrx(nq,nq), EVec(nq,nq), EVal(nq*(nq+1)/2), uMtrx(nX), Scrt(nq,nX)
logical g12K, Diagonal
real*8, parameter :: Zero_Approx = 0.1D-9
real*8, allocatable :: Work(:), W(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *

do i=1,nq
  do j=1,nx
    if (abs(Bmtrx(i,j)) < 1.0D-10) Bmtrx(i,j) = Zero
  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('ElRed: The B matrix','(5e21.12)',Bmtrx,nq,nx)
call RecPrt('ElRed: The u matrix','(5e21.12)',umtrx,nx,1)
#endif
if (nq == 0) then
  nK = 0
  Go To 99
end if

!                         T
! Form the G matrix, G=BuB

do j=1,nX
  do i=1,nq
    Scrt(i,j) = BMtrx(i,j)*umtrx(j)
  end do
end do
call DGEMM_('N','T',nq,nq,nX,1.0d0,Scrt,nq,Bmtrx,nq,0.0d0,Gmtrx,nq)

Diagonal = .true.
do i=1,nq
  Sum = 0.0d0
  do j=1,nq
    if (abs(Gmtrx(i,j)) < 1.0D-10) Gmtrx(i,j) = Zero
    if (j /= i) Sum = Sum+GMtrx(i,j)
  end do
  Diagonal = Diagonal .and. (Sum == 0.0d0)
end do

#ifdef _DEBUGPRINT_
call RecPrt('ElRed: The G Matrix (nq x nq)','(5e21.12)',Gmtrx,nq,nq)
write(6,*) 'Diagonal=',Diagonal
#endif

! Set up a unit matrix

call dcopy_(nq*nq,[Zero],0,EVec,1)
call dcopy_(nq,[One],0,EVec,nq+1)

! Set up the Hessian in lower triangular form, the elements are symmetrized.

do i=1,nQ
  do j=1,i
    ijTri = i*(i-1)/2+j
    EVal(ijTri) = Half*(Gmtrx(i,j)+Gmtrx(j,i))
  end do
end do
#ifdef _DEBUGPRINT_
call TriPrt('Eval prediagonalization',' ',EVal,nQ)
#endif

!                                                    |g 0|
! Compute eigenvalues and eigenvectors G(K L) = (K L)|0 0|
! K: nonredundant vectors with eigenvalues g
! L: redundant vectors

if (.not. Diagonal) then
  N = nQ
  LDZ = max(1,N)
  call mma_allocate(Work,3*N,Label='Work')
  Work(:) = Zero
  call mma_allocate(W,N,Label='W')
  W(:) = Zero
  Info = 0
  call dspev_('V','U',N,Eval,W,EVec,LDZ,Work,Info)
  if (Info /= 0) then
    write(6,*) 'Info /= 0'
    write(6,*) 'Info=',Info
    call Abend()
  end if
  call FZero(EVal,N*(N+1)/2)
  do i=1,N
    ii = i*(i+1)/2
    EVal(ii) = W(i)
  end do
  call mma_deallocate(W)
  call mma_deallocate(Work)
end if
call DScal_(nQ*(nQ+1)/2,-1.0d0,EVal,1)
call JacOrd(EVal,EVec,nQ,nQ)
!Fix standard direction.
do iQ=1,nQ
  call VecPhase(EVec(1,iQ),nQ)
end do
call DScal_(nQ*(nQ+1)/2,-1.0d0,EVal,1)
#ifdef _DEBUGPRINT_
call RecPrt('ElRed: Eigenvectors',' ',EVec,nQ,nQ)
call TriPrt('ElRed: Eigenvalues',' ',EVal,nQ)
#endif

!                                    -1/2
! Remove redundant vectors and form g     K

nK = 0
do i=1,nQ
  ii = i*(i+1)/2
  if (EVal(ii) > Thr) then
    nK = nK+1
  end if
  EVal(i) = EVal(ii)
  !if (g12K .and. (abs(EVal(i)) > Zero))
  if (g12K .and. (abs(EVal(i)) > Zero_Approx)) call DScal_(nQ,One/sqrt(EVal(i)),EVec(1,i),1)
end do
#ifdef _DEBUGPRINT_
call RecPrt('ElRed: The NonRedundant eigenvectors','(5e21.12)',EVec,nQ,nK)
call RecPrt('ElRed: eigenvalues ','(8E12.4)',EVal,1,nK)
#endif

99 continue

return

end subroutine ElRed
