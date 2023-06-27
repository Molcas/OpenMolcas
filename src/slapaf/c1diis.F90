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

use Slapaf_Parameters, only: iOptC

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8 q(nInter,nIter+1), dq(nInter,nIter), H(nInter,nInter), g(nInter,nIter+1), error(nInter,nIter), B((nIter+1)*(nIter+1)), &
       RHS(nIter+1)
integer iP(nIter), iRc
real*8, allocatable :: A(:)
! Statement function
ij(i,j,lda) = (j-1)*lda+i

iRout = 114
iPrint = nPrint(iRout)

call mma_allocate(A,nInter**2,Label='A')
call dcopy_(nInter**2,H,1,A,1)
iRc = 0
call dpotrf_('U',nInter,A,nInter,iRC)
if (iRC /= 0) then
  write(6,*) 'C1DIIS(DPOTRF): iRC=',iRC
  call Abend()
end if

! Compute the new set of error vectors

call dcopy_(nInter*nIter,g,1,Error,1)
iRc = 0
call DPOTRS('U',nInter,nIter,A,nInter,Error,nInter,iRC)
if (iRC /= 0) then
  write(6,*) 'C1DIIS(DPOTRS): iRC=',iRC
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
  if (iand(iOptC,16) == 16) then
    Err1 = DDot_(nInter,Error(1,iP(i)),1,Error(1,iP(i)),1)
  else if (iand(iOptC,32) == 32) then
    Err1 = DDot_(nInter,Error(1,iP(i)),1,g(1,iP(i)),1)
  else if (iand(iOptC,64) == 64) then
    Err1 = DDot_(nInter,g(1,iP(i)),1,g(1,iP(i)),1)
  else
    Err1 = Zero
    call WarningMessage(2,' Illegal iOptC setting!')
    call Abend()
  end if
  ii = i
  do j=i+1,nIter
    if (iand(iOptC,16) == 16) then
      Err2 = DDot_(nInter,Error(1,iP(j)),1,Error(1,iP(j)),1)
    else if (iand(iOptC,32) == 32) then
      Err2 = DDot_(nInter,Error(1,iP(j)),1,g(1,iP(j)),1)
    else if (iand(iOptC,64) == 64) then
      Err2 = DDot_(nInter,g(1,iP(j)),1,g(1,iP(j)),1)
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
if (iPrint >= 99) write(6,*) ' iP=',iP

MaxWdw = max(2,(nInter-nFix)/2)
mIter = min(nIter,min(MinWdw,MaxWdw))
iOff = max(0,nIter-mIter)
B(ij(mIter+1,mIter+1,mIter+1)) = Zero
RHS(mIter+1) = -One
do i=1,mIter
  do j=1,i-1
    if (iand(iOptC,16) == 16) then
      B(ij(i,j,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,error(1,iP(j+iOff)),1)
    else if (iand(iOptC,32) == 32) then
      B(ij(i,j,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,g(1,iP(j+iOff)),1)
    else if (iand(iOptC,64) == 64) then
      B(ij(i,j,mIter+1)) = DDot_(nInter,g(1,iP(i+iOff)),1,g(1,iP(j+iOff)),1)
    else
      call WarningMessage(2,' Illegal iOptC setting!')
      call Abend()
    end if
    B(ij(j,i,mIter+1)) = B(ij(i,j,mIter+1))
  end do
  if (iand(iOptC,16) == 16) then
    B(ij(i,i,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,error(1,iP(i+iOff)),1)
  else if (iand(iOptC,32) == 32) then
    B(ij(i,i,mIter+1)) = DDot_(nInter,error(1,iP(i+iOff)),1,g(1,iP(i+iOff)),1)
  else if (iand(iOptC,64) == 64) then
    B(ij(i,i,mIter+1)) = DDot_(nInter,g(1,iP(i+iOff)),1,g(1,iP(i+iOff)),1)
  else
    call WarningMessage(2,' Illegal iOptC setting!')
    call Abend()
  end if
  B(ij(i,mIter+1,mIter+1)) = -One
  B(ij(mIter+1,i,mIter+1)) = -One
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

call dcopy_(nInter,[Zero],0,q(1,nIter+1),1)
call dcopy_(nInter,[Zero],0,g(1,nIter+1),1)
do jIter=1,mIter
  call DaXpY_(nInter,RHS(jIter),q(1,iP(jIter+iOff)),1,q(1,nIter+1),1)
  call DaXpY_(nInter,RHS(jIter),g(1,iP(jIter+iOff)),1,g(1,nIter+1),1)
end do
if (iPrint >= 99) then
  call RecPrt(' The ipv',' ',q(1,nIter+1),1,nInter)
  call RecPrt(' The igv',' ',g(1,nIter+1),1,nInter)
end if

! Compute a new independent geometry by relaxation of
! the interpolated gradient vector.

call dcopy_(nInter,g(1,nIter+1),1,dq(1,nIter),1)
call DPOTRS('U',nInter,1,A,nInter,dq(1,nIter),nInter,iRC)
if (iRC /= 0) then
  write(6,*) 'C1DIIS(DPOTRS): iRC=',iRC
  call Abend()
end if
if (iPrint >= 99) call RecPrt(' dq',' ',dq(1,nIter),1,nInter)

! The shift is relative to the interpolated parameter
! vector and we have to change it so that it is relative to the
! actual parameter vector.

do iInter=1,nInter
  dq(iInter,nIter) = dq(iInter,nIter)+q(iInter,nIter+1)-q(iInter,nIter)
end do
if (iPrint >= 99) call RecPrt(' dq(corr.)',' ',dq(1,nIter),1,nInter)

call mma_deallocate(A)

return

end subroutine C1DIIS
