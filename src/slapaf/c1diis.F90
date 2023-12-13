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
! Copyright (C) 1994, Roland Lindh                                     *
!***********************************************************************

subroutine C1DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,nFix,iP,MinWdw)
!***********************************************************************
!                                                                      *
!         References:                                                  *
!           C1-DIIS: P. Csaszar and P. Pulay, J. Mol. Struc.           *
!                    114, 31-34 (1984).                                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December 1994                                            *
!***********************************************************************

use Slapaf_Info, only: iOptC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nInter, nIter, nFix, MinWdw
real(kind=wp), intent(inout) :: q(nInter,nIter+1), dq(nInter,nIter), g(nInter,nIter+1)
real(kind=wp), intent(in) :: H(nInter,nInter)
real(kind=wp), intent(out) :: error(nInter,nIter), B((nIter+1)**2), RHS(nIter+1)
integer(kind=iwp), intent(out) :: iP(nIter)
#include "print.fh"
integer(kind=iwp) :: i, ii, ij, iOff, iPrint, iRc, iRout, iSave, j, jIter, MaxWdw, mIter
real(kind=wp) :: Err1, Err2
real(kind=wp), allocatable :: A(:,:)
real(kind=wp), external :: DDot_

iRout = 114
iPrint = nPrint(iRout)

call mma_allocate(A,nInter,nInter,Label='A')
A(:,:) = H(:,:)
iRc = 0
call dpotrf_('U',nInter,A,nInter,iRC)
if (iRC /= 0) then
  write(u6,*) 'C1DIIS(DPOTRF): iRC=',iRC
  call Abend()
end if

! Compute the new set of error vectors

Error(:,:) = g(:,1:nIter)
iRc = 0
call DPOTRS('U',nInter,nIter,A,nInter,Error,nInter,iRC)
if (iRC /= 0) then
  write(u6,*) 'C1DIIS(DPOTRS): iRC=',iRC
  call Abend()
end if
if (iPrint >= 99) call RecPrt(' Error vectors',' ',error,nInter,nIter)

! Set up small system of linear equations
! If more error vectors than degrees of freedom
! exclude those with large error.

do i=1,nIter
  iP(i) = i
end do

! Bubble sort index array with respect to the magnitude of the
! error vector.

do i=1,nIter-1
  if (btest(iOptC,4)) then
    Err1 = DDot_(nInter,Error(:,iP(i)),1,Error(:,iP(i)),1)
  else if (btest(iOptC,5)) then
    Err1 = DDot_(nInter,Error(:,iP(i)),1,g(:,iP(i)),1)
  else if (btest(iOptC,6)) then
    Err1 = DDot_(nInter,g(:,iP(i)),1,g(:,iP(i)),1)
  else
    Err1 = Zero
    call WarningMessage(2,' Illegal iOptC setting!')
    call Abend()
  end if
  ii = i
  do j=i+1,nIter
    if (btest(iOptC,4)) then
      Err2 = DDot_(nInter,Error(:,iP(j)),1,Error(:,iP(j)),1)
    else if (btest(iOptC,5)) then
      Err2 = DDot_(nInter,Error(:,iP(j)),1,g(:,iP(j)),1)
    else if (btest(iOptC,6)) then
      Err2 = DDot_(nInter,g(:,iP(j)),1,g(:,iP(j)),1)
    else
      Err2 = Zero
      call WarningMessage(2,' Illegal iOptC setting!')
      call Abend()
    end if
    if (Err2 > Err1) then
      ii = j
      Err1 = Err2
    end if
  end do
  if (ii /= i) then
    iSave = iP(i)
    iP(i) = iP(ii)
    iP(ii) = iSave
  end if
end do
if (iPrint >= 99) write(u6,*) ' iP=',iP

MaxWdw = max(2,(nInter-nFix)/2)
mIter = min(nIter,min(MinWdw,MaxWdw))
iOff = max(0,nIter-mIter)
B((mIter+1)**2) = Zero
RHS(mIter+1) = -One
do i=1,mIter
  do j=1,i
    ij = (j-1)*(mIter+1)+i
    if (btest(iOptC,4)) then
      B(ij) = DDot_(nInter,error(:,iP(i+iOff)),1,error(:,iP(j+iOff)),1)
    else if (btest(iOptC,5)) then
      B(ij) = DDot_(nInter,error(:,iP(i+iOff)),1,g(:,iP(j+iOff)),1)
    else if (btest(iOptC,6)) then
      B(ij) = DDot_(nInter,g(:,iP(i+iOff)),1,g(:,iP(j+iOff)),1)
    else
      call WarningMessage(2,' Illegal iOptC setting!')
      call Abend()
    end if
    if (j < i) B((i-1)*(mIter+1)+j) = B(ij)
  end do
  B(mIter*(mIter+1)+i) = -One
  B(i*(mIter+1)) = -One
  RHS(i) = Zero
end do
if (iPrint >= 99) then
  call RecPrt(' The B Matrix',' ',B,mIter+1,mIter+1)
  call RecPrt(' The RHS',' ',RHS,1,mIter+1)
end if

! Solve linear equation system

call Gauss(mIter+1,mIter+1,B,RHS,RHS)
if (iPrint >= 99) call RecPrt(' The solution vector',' ',RHS,1,mIter+1)

! Compute the interpolated parameter vector and
! the interpolated gradient vector.

q(:,nIter+1) = Zero
g(:,nIter+1) = Zero
do jIter=1,mIter
  q(:,nIter+1) = q(:,nIter+1)+RHS(jIter)*q(:,iP(jIter+iOff))
  g(:,nIter+1) = g(:,nIter+1)+RHS(jIter)*g(:,iP(jIter+iOff))
end do
if (iPrint >= 99) then
  call RecPrt(' The ipv',' ',q(:,nIter+1),1,nInter)
  call RecPrt(' The igv',' ',g(:,nIter+1),1,nInter)
end if

! Compute a new independent geometry by relaxation of
! the interpolated gradient vector.

dq(:,nIter) = g(:,nIter+1)
call DPOTRS('U',nInter,1,A,nInter,dq(:,nIter),nInter,iRC)
if (iRC /= 0) then
  write(u6,*) 'C1DIIS(DPOTRS): iRC=',iRC
  call Abend()
end if
if (iPrint >= 99) call RecPrt(' dq',' ',dq(:,nIter),1,nInter)

! The shift is relative to the interpolated parameter
! vector and we have to change it so that it is relative to the
! actual parameter vector.

dq(:,nIter) = dq(:,nIter)+q(:,nIter+1)-q(:,nIter)
if (iPrint >= 99) call RecPrt(' dq(corr.)',' ',dq(:,nIter),1,nInter)

call mma_deallocate(A)

return

end subroutine C1DIIS
