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

subroutine Eq_Solver(Mode,M,N,NRHS,B,Curvilinear,Degen,dSS,DFC)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
character, intent(in) :: Mode
integer(kind=iwp), intent(in) :: M, N, NRHS
real(kind=wp), intent(in) :: B(M,N), Degen(M), dSS(*)
logical(kind=iwp), intent(in) :: Curvilinear
real(kind=wp), intent(_OUT_) :: DFC(*)
#include "warnings.h"
integer(kind=iwp) :: i, ij, INFO, iRHS, LDA, LDB, LWork
real(kind=wp) :: Temp(1)
real(kind=wp), allocatable :: A(:,:), Btmp(:,:), Work(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Solve the equation Ax=b

LDA = M
LDB = max(1,M,N)
if (Mode == 'T') then
  call mma_allocate(A,M,M,Label='A')
  A(:,1:N) = B(:,:)
  A(:,N+1:) = Zero
  if (.not. Curvilinear) then
    do i=1,M
      A(i,:) = sqrt(Degen(i))*A(i,:)
    end do
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('A',' ',A,M,M)
# endif
else
  call mma_allocate(A,M,N,Label='A')
  A(:,:) = B(:,:)
  if (.not. Curvilinear) then
    do i=1,M
      A(i,:) = sqrt(Degen(i))*A(i,:)
    end do
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('A',' ',A,M,N)
# endif
end if

call mma_allocate(Btmp,LDB,NRHS,Label='Btmp')
Btmp(:,:) = Zero

if (Mode == 'T') then
  do iRHS=1,nRHS
    ij = (iRHS-1)*N
    Btmp(1:N,iRHS) = dss(ij+1:ij+N)
  end do
else
# ifdef _DEBUGPRINT_
  call RecPrt('B(raw)',' ',dss,M,nRHS)
# endif
  do iRHS=1,nRHS
    if (.not. Curvilinear) then
      ij = (iRHS-1)*M
      Btmp(1:M,iRHS) = dss(ij+1:ij+M)*sqrt(Degen(:))
    else
      ij = (iRHS-1)*M
      Btmp(1:M,iRHS) = dss(ij+1:ij+M)
    end if
  end do
end if
#ifdef _DEBUGPRINT_
call RecPrt('B(in)',' ',Btmp,LDB,NRHS)
#endif

LWork = -1
INFO = 0
call dgels_(Mode,M,N,NRHS,A,LDA,Btmp,LDB,Temp,LWork,INFO)
LWork = int(Temp(1))
call mma_allocate(Work,LWork,Label='Work')
INFO = 0
call dgels_(Mode,M,N,NRHS,A,LDA,Btmp,LDB,Work,LWork,INFO)
if (INFO > 0) then
  call WarningMessage(2,'Error in Eq_Solver')
  write(u6,*)
  write(u6,*) '***********************************************'
  write(u6,*) ' ERROR: Eq_Solver could not find a solution.   '
  write(u6,*) ' The matrix is rank deficient.                 '
  write(u6,*) '***********************************************'
  call Quit(_RC_INTERNAL_ERROR_)
end if

#ifdef _DEBUGPRINT_
call RecPrt('B(out)',' ',Btmp,LDB,NRHS)
#endif
if (Mode == 'T') then
  do iRHS=1,nRHS
    if (.not. Curvilinear) Btmp(1:M,iRHS) = Btmp(1:M,iRHS)/sqrt(Degen(:))
    ij = (iRHS-1)*M
    DFC(ij+1:ij+M) = Btmp(1:M,iRHS)
  end do
else
  do iRHS=1,nRHS
    ij = (iRHS-1)*N
    DFC(ij+1:ij+N) = Btmp(1:N,iRHS)
  end do
end if

call mma_deallocate(Work)
call mma_deallocate(Btmp)
call mma_deallocate(A)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Eq_Solver
