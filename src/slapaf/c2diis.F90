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
! Copyright (C) 1994,1995, Roland Lindh                                *
!***********************************************************************

subroutine C2DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,Scrt1,nScrt1,nFix,iP)
!***********************************************************************
!                                                                      *
!         References:                                                  *
!           C2-DIIS: Sellers, Int. J. Quantum Chem. 45, 31-41(1993).   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             December '94                                             *
!                                                                      *
!             Modified for anharmonic constants by R. Lindh, Oct. '95  *
!***********************************************************************

use Slapaf_Parameters, only: iOptC

implicit real*8(a-h,o-z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8 q(nInter,nIter+1), dq(nInter,nIter), H(nInter,nInter), g(nInter,nIter+1), error(nInter,nIter+1), B((nIter+1)*(nIter+1)), &
       RHS(nIter+1), Scrt1(nScrt1)
integer iP(nIter), iRc
logical Fail
real*8, allocatable :: A(:)
! Statement function
ij(i,j,lda) = (j-1)*lda+i
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

iRout = 121
iPrint = nPrint(iRout)

call FZero(Error,nInter*(nIter+1))

call mma_allocate(A,nInter**2,Label='A')
call dcopy_(nInter**2,H,1,A,1)
iRc = 0
call dpotrf_('U',nInter,A,nInter,iRC)
if (iRC /= 0) then
  write(6,*) 'C2DIIS(DPOTRF): iRC=',iRC
  call Abend()
end if

! Compute the new set of error vectors. Here we have two
! possibilities depending on if we use a 2nd or 3rd order
! update method.
!
! Note!!!!
!
! i)  We store the force in g Not the gradient
!
! ii) e is the displacement which should be added to the
!     current position to get to the equilibrium geometry.
!
!    r  = r   + e      e   =r  - r
!     eq   i-1   i-1    i-1  eq   i-1
!
! Compute:
!    r  = r + e        e =r  - r
!     eq   i   i        k  eq   k

if ((iand(iOptC,16) == 16) .or. (iand(iOptC,32) == 32)) then
  call ThrdO(nInter,g(1,nIter),A,Error(1,nIter),Fail)
  if (Fail) then
    call WarningMessage(2,'C2Diis: ThrdO Failed!')
    call Quit_OnConvError()
  end if
  do iIter=1,nIter-1
    do iInter=1,nInter
      Error(iInter,iIter) = Error(iInter,nIter)+q(iInter,nIter)-q(iInter,iIter)
    end do
  end do
end if
if (iPrint >= 99) call RecPrt(' Error vectors',' ',error,nInter,nIter)

! Set up small system of linear equations
! If more error vectors than degrees of freedom
! exclude those with large error.

do i=1,nIter
  iP(i) = i
end do

! Bubble sort index array with respect to the magnitude of
!
! <g|g>
! <g|dx>
! <dx|g>
! <dx|dx>

do i=max(1,nIter-11),nIter-1
  if (iand(iOptC,16) == 16) then
    Err1 = DDot_(nInter,Error(1,iP(i)),1,Error(1,iP(i)),1)
  else if (iand(iOptC,32) == 32) then
    Err1 = DDot_(nInter,g(1,iP(i)),1,Error(1,iP(i)),1)
  else if (iand(iOptC,64) == 64) then
    Err1 = DDot_(nInter,g(1,iP(i)),1,g(1,iP(i)),1)
  else
    Err1 = Zero
    call WarningMessage(2,' Illegal iOptC setting!')
    call Quit_OnUserError()
  end if
  ii = i
  do j=i+1,nIter
    if (iand(iOptC,16) == 16) then
      Err2 = DDot_(nInter,Error(1,iP(j)),1,Error(1,iP(j)),1)
    else if (iand(iOptC,32) == 32) then
      Err2 = DDot_(nInter,g(1,iP(j)),1,Error(1,iP(j)),1)
    else if (iand(iOptC,64) == 64) then
      Err2 = DDot_(nInter,g(1,iP(j)),1,g(1,iP(j)),1)
    else
      Err2 = Zero
      call WarningMessage(2,' Illegal iOptC setting!')
      call Quit_OnUserError()
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

!MaxWdw = max(3,3*(nInter-nFix)/4)
MaxWdw = max(3,(nInter-nFix)/2)
MinWdw = min(5,MaxWdw)
mIter = min(nIter,MinWdw)
iOff = max(0,nIter-mIter)
Thrhld = 0.1D-13
ThrCff = (Two*Ten)**2
ThrLdp = Ten**3
do i=1,mIter
  do j=1,i
    if (iand(iOptC,16) == 16) then
      B(iTri(i,j)) = DDot_(nInter,Error(1,iP(i+iOff)),1,Error(1,iP(j+iOff)),1)
    else if (iand(iOptC,32) == 32) then
      B(iTri(i,j)) = DDot_(nInter,Error(1,iP(i+iOff)),1,g(1,iP(j+iOff)),1)
    else if (iand(iOptC,64) == 64) then
      B(iTri(i,j)) = DDot_(nInter,g(1,iP(i+iOff)),1,g(1,iP(j+iOff)),1)
    else
      call WarningMessage(2,' Illegal iOptC setting!')
      call Quit_OnUserError()
    end if
  end do
end do
if (iPrint >= 99) call TriPrt(' The B Matrix',' ',B,mIter)
call dcopy_(mIter**2,[Zero],0,Scrt1,1)
call dcopy_(mIter,[One],0,Scrt1,mIter+1)
call NIDiag_new(B,Scrt1,mIter,mIter)
if (iPrint >= 99) then
  call TriPrt(' The B Matrix after diagonalization','(9E10.2)',B,mIter)
  call RecPrt(' Eigenvectors','(9E10.2)',Scrt1,mIter,mIter)
end if

! Renormalize the eigenvectors and eigenvalues to the
! C1-DIIS format.

do iVec=1,mIter
  Alpha = Zero
  do i=1,mIter
    Alpha = Alpha+Scrt1(ij(i,iVec,mIter))
  end do
  Alpha = One/Alpha
  call DScal_(mIter,Alpha,Scrt1(ij(1,iVec,mIter)),1)
  B(iTri(iVec,iVec)) = B(iTri(iVec,iVec))*Alpha**2
end do
if (iPrint >= 99) then
  write(6,*) ' After normalization to C1-DIIS format'
  call TriPrt(' The B Matrix after diagonalization','(9E10.2)',B,mIter)
  call RecPrt(' Eigenvectors',' ',Scrt1,mIter,mIter)
end if

! Select a vector.

ee_old = 1.0D+72
c2_old = 1.0D+72
iVec_old = -99999999
do iVec=1,mIter
  if (iPrint >= 99) write(6,*) ' Scanning vector',iVec
  ee_new = B(iTri(iVec,iVec))
  if (iPrint >= 99) write(6,*) ' ee_old, ee_new=',ee_old,ee_new

  ! Examine if <e|e> is too low (possible round-off) or linear dependency.

  if (ee_new < Thrhld) then
    if (iPrint >= 99) write(6,*) ' <e|e> is low in DIIS, iVec,<e|e>=',iVec,ee_new

    ! Reject if coefficients are too large (linear dep.).

    c2_new = DDot_(mIter,Scrt1(ij(1,iVec,mIter)),1,Scrt1(ij(1,iVec,mIter)),1)
    if (c2_new > ThrCff) then
      if (iPrint >= 99) write(6,*) ' c**2 is too large in DIIS, iVec,c**2=',iVec,c2_new
      Go To 10
    end if
  end if

  ! Reject if coefficients are by far too large (linear dep.).

  c2_new = DDot_(mIter,Scrt1(ij(1,iVec,mIter)),1,Scrt1(ij(1,iVec,mIter)),1)
  if (c2_new > ThrLdp) then
    if (iPrint >= 99) write(6,*) ' c**2 is too large in DIIS, iVec,c**2=',iVec,c2_new
    Go To 10
  end if

  ! Keep the best candidate

  if (ee_new*Five < ee_old) then
    ! New vector much lower eigenvalue.
    c2_old = c2_new
    ee_old = ee_new
    iVec_old = iVec
    if (iPrint >= 99) write(6,*) 'New vector much lower eigenvalue',iVec_old
  else if (ee_new <= ee_old*Five) then
    ! New vector is close to the old vector.
    ! Selection based on relative weight of the last geometry.
    if (iPrint >= 99) write(6,*) 'Eigenvalues are close',iVec_old,iVec
    t1 = abs(Scrt1(ij(mIter,iVec_old,mIter)))/sqrt(c2_old)
    t2 = abs(Scrt1(ij(mIter,iVec,mIter)))/sqrt(c2_new)
    if (t2 > t1*1.2d0) then
      ! New vector much better relative weight.
      c2_old = c2_new
      ee_old = ee_new
      iVec_old = iVec
      if (iPrint >= 99) write(6,*) 'New vector much better relative weight',iVec_old
    else if (t2*1.2d0 < t1) then
      ! Vectors are close in relative weight too!
      ! Select on eigenvalue only
      if (iPrint >= 99) write(6,*) 'Relative weights are close',iVec_old,iVec
      if (ee_new < ee_old) then
        c2_old = c2_new
        ee_old = ee_new
        iVec_old = iVec
        if (iPrint >= 99) write(6,*) 'New vector has lower eigenvalue',iVec_old
      end if
    end if
  end if

10 continue
end do
if ((iVec_old < 1) .or. (iVec_old > mIter)) then
  call WarningMessage(2,' No proper solution found in C2-DIIS!')
  call Abend()
end if
call dcopy_(mIter,Scrt1(ij(1,iVec_old,mIter)),1,RHS,1)

if (iPrint >= 99) then
  write(6,*) ' Selecting root',iVec_old
  call RecPrt(' The solution vector',' ',RHS,1,mIter)
end if

! Compute the interpolated parameter vector and
! the interpolated gradient vector.

call dcopy_(nInter,[Zero],0,q(1,nIter+1),1)
call dcopy_(nInter,[Zero],0,g(1,nIter+1),1)
call dcopy_(nInter,[Zero],0,Scrt1,1)
do iIter=1,mIter

  if (abs(RHS(iIter)) < 1.0D-12) Go To 11

  call DaXpY_(nInter,RHS(iIter),Error(1,iP(iIter+iOff)),1,Scrt1,1)

  ! The interpolated parameter vector is computed as a
  ! simple linear combination

  call DaXpY_(nInter,RHS(iIter),q(1,iP(iIter+iOff)),1,q(1,nIter+1),1)

  ! The interpolated gradient vector (Stored as force)
  !
  ! Sum(i) c  g                 (Coeffs stored in RHS)
  !         i  i

  call DaXpY_(nInter,RHS(iIter),g(1,iP(iIter+iOff)),1,g(1,nIter+1),1)

11 continue
end do

if (iPrint >= 99) then
  call RecPrt(' The iev',' ',Scrt1,1,nInter)
  call RecPrt(' The ipv',' ',q(1,nIter+1),1,nInter)
  call RecPrt(' The igv',' ',g(1,nIter+1),1,nInter)
end if

! Compute a new independent geometry by relaxation of
! the interpolated gradient vector.

call dcopy_(nInter,g(1,nIter+1),1,dq(1,nIter),1)
iRc = 0
call DPOTRS('U',nInter,1,A,nInter,dq(1,nIter),nInter,iRC)
if (iRC /= 0) then
  write(6,*) 'C2DIIS(DPOTRS): iRC=',iRC
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

end subroutine C2DIIS
