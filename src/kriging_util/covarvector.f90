
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
SUBROUTINE covarVector(gh)
  use kriging_mod
  Implicit None
#include "stdalloc.fh"
  integer i,i0,i1,j,j0,j1,k,k0,k1,gh
  real*8 sdiffxi,sdiffxj,sdiffxk
  real*8, Allocatable ::  diffxi(:),diffxj(:), diffxk(:)
!#define _DEBUGPRINT_

  Call mma_Allocate(diffxi,nPoints,label="diffxi")
  Call mma_Allocate(diffxj,nPoints,label="diffxj")
  Call mma_Allocate(diffxk,nPoints,label="diffxk")
!
  cv = 0
  i0 = 0
  call defdlrl()
!
! Covariant Vector in kriging - First part of eq (4) in ref.
!
  if (gh.eq.0) then
!
    call matern(dl, cv(1:nPoints,1,1), nPoints,1)
    call matderiv(1, dl, cvMatFDer, nPoints, 1)
    do i=1,nInter
!     1st derivatives second part of eq. (4)
      diffxi(:) = 2.0D0*rl(:,i)/l(i)
      i0 = nPoints + 1 + (i-1)*(nPoints-nD)
      i1 = i0 + (nPoints-nD) - 1
      cv(i0:i1,1,1) = cvMatFder(1+nD:nPoints) * diffxi(1+nD:nPoints)
    enddo
#ifdef _DEBUGPRINT_
    Call RecPrt(' The covector for energies','(12(2x,E9.3))',cv(:,1,1),m_t,1)
#endif
! Covariant vector in Gradient Enhanced Kriging
!
  else if(gh.eq.1) then
!
    call matderiv(1, dl, cvMatFder, nPoints, 1)
    call matderiv(2, dl, cvMatSder, nPoints, 1)
    do i=1,nInter
      diffxi(:) = 2.0D0*rl(:,i)/l(i)
      cv(1:nPoints,i,1) = -cvMatFder(1:nPoints) * diffxi(1:nPoints)
      do j = 1,nInter
        j0 = nPoints + 1 + (j-1)*(nPoints-nD)
        j1 = j0 + (nPoints-nD) - 1
        diffxj(:) = -2.0D0*rl(:,j)/l(j)
        if (i.eq.j) Then
         cv(j0:j1,i,1) = cvMatSder(1+nD:nPoints) * diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints) &
                       - cvMatFder(1+nD:nPoints)*(2/(l(i)*l(j)))
        else
         cv(j0:j1,i,1) = cvMatSder(1+nD:nPoints) * diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)
        end if
      enddo
    enddo
#ifdef _DEBUGPRINT_
    Call RecPrt(' The covector for gradients','(12(2x,E9.3))',cv(:,:,1),m_t,nInter)
#endif
!
  else if(gh.eq.2) then
!
      !    print *,'covar vector calling deriv(3) for Kriging Hessian'
    call matderiv(1, dl, cvMatFder, nPoints, 1)
    call matderiv(2, dl, cvMatSder, nPoints, 1)
    call matderiv(3, dl, cvMatTder, nPoints, 1)
    do i = 1, nInter
      diffxi(:) = 2.0D0*rl(:,i)/l(i)
      sdiffxi = 2.0D0/l(i)**2
      do j = 1, nInter
        diffxj(:) = 2.0D0*rl(:,j)/l(j)
        sdiffxj = 2.0D0/l(j)**2
        if (i.eq.j) Then
          cv(1:nPoints,i,j) = cvMatSder(:) * diffxi(:)*diffxj(:) + cvMatFder(:)*2.0D0/(l(i)*l(j))
        else
          cv(1:nPoints,i,j) = cvMatSder(:) * diffxi(:)*diffxj(:)
        end if
        do k = 1, nInter
          diffxk(:) = 2.0D0*rl(:,k)/l(k)
          sdiffxk = 2.0D0/l(k)**2
          k0 = nPoints + 1 + (k-1)*(nPoints-nD)
          k1 = k0 + (nPoints-nD) - 1
          if (i.eq.j.and.j.eq.k) then
            cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints) &
                    + 3.0D0*cvMatSder(1+nD:nPoints)*diffxi(1+nD:nPoints)*sdiffxj
          else if (i.eq.j) then
            cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints) &
                          + cvMatSder(1+nD:nPoints)*diffxk(1+nD:nPoints)*sdiffxi
          else if (i.eq.k) then
            cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints) &
                          + cvMatSder(1+nD:nPoints)*diffxj(1+nD:nPoints)*sdiffxi
          else if (j.eq.k) then
            cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints) &
                          + cvMatSder(1+nD:nPoints)*diffxi(1+nD:nPoints)*sdiffxk
          else
            cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints)
          endif
        enddo
      enddo
    enddo
  else
    Write (6,*) ' Illegal value of gh:',gh
    Call Abend()
  endif
!
  Call mma_deallocate(diffxi)
  Call mma_deallocate(diffxj)
  Call mma_deallocate(diffxk)
!
contains
!
SUBROUTINE defdlrl()
  use kriging_mod
  integer i,j

  dl(:)=0.0D0
  do i=1,nInter
    do j=1,nPoints
       rl(j,i) = (x(i,j) - x0(i))/l(i)
    enddo
    dl(:) = dl(:) + rl(:,i)**2
  enddo
END Subroutine defdlrl
!
END Subroutine covarvector
