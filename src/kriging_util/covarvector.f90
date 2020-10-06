
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
  integer i,i0,i1,j,j0,j1,k,k0,k1,gh,nInter
  real*8 sdiffxi,sdiffxj,sdiffxk
  real*8, Allocatable ::  diffxi(:),diffxj(:), diffxk(:)

  nInter =nInter_save

  Call mma_Allocate(diffxi,nPoints_v,label="diffxi")
  Call mma_Allocate(diffxj,nPoints_v,label="diffxj")
  Call mma_Allocate(diffxk,nPoints_v,label="diffxk")
!
  cv = 0
  i0 = 0
  call defdlrl()
!
! Covariant Vector in kriging - First part of eq (4) in ref.
!
  if (gh.eq.0) then
!
    call matern(dl, cv(1:nPoints_v,1,1,1), nPoints_v,1)
    call matderiv(1, dl, cvMatFDer, nPoints_v, 1)
    do i=1,nInter
!     1st derivatives second part of eq. (4)
      diffxi(:) = 2.0D0*rl(:,i)/l(i)
      i0 = nPoints_v + 1 + (i-1)*nPoints_g
      i1 = i0 + nPoints_g - 1
      cv(i0:i1,1,1,1) = cvMatFder(:,1) * diffxi(:)
    enddo
! Covariant vector in Gradient Enhanced Kriging
!
  else if(gh.eq.1) then
!
    call matderiv(1, dl, cvMatFder, nPoints_v, 1)
    call matderiv(2, dl, cvMatSder, nPoints_v, 1)
    do i=1,nInter
      diffxi(:) = 2.0D0*rl(:,i)/l(i)
      cv(1:nPoints_v,1,i,1) = -cvMatFder(:,1) * diffxi(:)
      do j = 1,nInter
        j0 = nPoints_v + 1 + (j-1)*nPoints_g
        j1 = j0 + nPoints_g - 1
        diffxj(:) = -2.0D0*rl(:,j)/l(j)
        if (i.eq.j) Then
         cv(j0:j1,1,i,1) = cvMatSder(:,1) * diffxi(:)*diffxj(:) - cvMatFder(:,1)*(2/(l(i)*l(j)))
        else
         cv(j0:j1,1,i,1) = cvMatSder(:,1) * diffxi(:)*diffxj(:)
        end if
      enddo
    enddo
!
  else if(gh.eq.2) then
!
      !    print *,'covar vector calling deriv(3) for Kriging Hessian'
    call matderiv(1, dl, cvMatFder, nPoints_v, 1)
    call matderiv(2, dl, cvMatSder, nPoints_v, 1)
    call matderiv(3, dl, cvMatTder, nPoints_v, 1)
    do i = 1, nInter
      diffxi(:) = 2.0D0*rl(:,i)/l(i)
      sdiffxi = 2.0D0/l(i)**2
      do j = 1, nInter
        diffxj(:) = 2.0D0*rl(:,j)/l(j)
        sdiffxj = 2.0D0/l(j)**2
        if (i.eq.j) Then
          cv(1:nPoints_v,1,i,j) = cvMatSder(:,1) * diffxi(:)*diffxj(:) + cvMatFder(:,1)*2.0D0/(l(i)*l(j))
        else
          cv(1:nPoints_v,1,i,j) = cvMatSder(:,1) * diffxi(:)*diffxj(:)
        end if
        do k = 1, nInter
          diffxk(:) = 2.0D0*rl(:,k)/l(k)
          sdiffxk = 2.0D0/l(k)**2
          k0 = nPoints_v + 1 + (k-1)*nPoints_g
          k1 = k0 + nPoints_g - 1
          if (i.eq.j.and.j.eq.k) then
            cv(k0:k1,1,i,j) = cvMatTder(:,1)*diffxi(:)*diffxj(:)*diffxk(:) + 3.0D0*cvMatSder(:,1)*diffxi(:)*sdiffxj
          else if (i.eq.j) then
            cv(k0:k1,1,i,j) = cvMatTder(:,1)*diffxi(:)*diffxj(:)*diffxk(:) + cvMatSder(:,1)*diffxk(:)*sdiffxi
          else if (i.eq.k) then
            cv(k0:k1,1,i,j) = cvMatTder(:,1)*diffxi(:)*diffxj(:)*diffxk(:) + cvMatSder(:,1)*diffxj(:)*sdiffxi
          else if (j.eq.k) then
            cv(k0:k1,1,i,j) = cvMatTder(:,1)*diffxi(:)*diffxj(:)*diffxk(:) + cvMatSder(:,1)*diffxi(:)*sdiffxk
          else
            cv(k0:k1,1,i,j) = cvMatTder(:,1)*diffxi(:)*diffxj(:)*diffxk(:)
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
  integer i,j,nInter

  nInter=nInter_save

  dl(:)=0.0D0
  do i=1,nInter
    do j=1,nPoints_v
       rl(j,i) = (x(i,j) - x0(i))/l(i)
    enddo
    dl(:) = dl(:) + rl(:,i)**2
  enddo
END Subroutine defdlrl
!
END Subroutine covarvector
