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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
SUBROUTINE covarMatrix(nPoints,nInter)
  use kriging_mod
  Implicit None
#include "stdalloc.fh"
  integer i,j,i0,i1,j0,j1,k,nPoints,nInter
  Real*8, Allocatable :: diffx_j(:,:), diffx_i(:,:), matFder(:,:),&
                         matSder(:,:), r(:,:,:), d(:,:)
!#define _DEBUG_
!
  Call mma_Allocate(diffx_j,nPoints,nPoints,Label="diffx_j")
  Call mma_Allocate(diffx_i,nPoints,nPoints,Label="diffx_i")
  Call mma_Allocate(matFder,nPoints,nPoints,Label="matFder")
  Call mma_Allocate(matSder,nPoints,nPoints,Label="matSder")
  Call mma_Allocate(r,nPoints,nPoints,nInter,Label="r")
  Call mma_Allocate(d,nPoints,nPoints,Label="d")
!
  full_R = 0
  d = 0
  diffx_j = 0
  diffx_i = 0
!
  do i=1,nInter
    do k=1,nPoints
      do j=1,nPoints
        r(k,j,i)=(x(i,k)-x(i,j))/l(i)
      end do
    end do
    d(:,:) = d(:,:) + r(:,:,i)**2
#ifdef _DEBUG_
    Call RecPrt('r',' ',r(1,1,i),nPoints,nPoints)
#endif
  end do
#ifdef _DEBUG_
  Call RecPrt('l',' ',l,1,nInter)
  Call RecPrt('x',' ',x,nInter,nPoints)
  Call RecPrt('d',' ',d,nPoints,nPoints)
#endif
!
    !Matern Function
  Call matern     (d, full_R(1:nPoints,1:nPoints), nPoints, nPoints)
!
! Writing the covariant matrix in GEK (eq 2 of DOI 10.1007/s00366-015-0397)
!
    !Matern first derivative with respect to d
  call matderiv(1, d, MatFder, nPoints, nPoints)
! Covariant matrix in Gradient Enhanced Kriging (eq 2 of DOI 10.1007/s00366-015-0397)):
!
    ! First line and first column derivative in Psi matrix
  do i=1,nInter
    i0 = i*nPoints +1
    i1 = i0        +nPoints-1
!
    diffx_i(:,:) = -2.0D0*r(:,:,i)/l(i)
!
    !  Writing the 1st row of 1st derivatives with respect the coordinates
    full_R(1:nPoints,i0:i1) = matFDer*diffx_i
    !  Writing the column of derivatives
    full_R(i0:i1,1:nPoints) = transpose(full_R(1:nPoints,i0:i1))
  enddo
!
    !Matern second derivative with respect to d
  call matderiv(2, d, matSder, nPoints, nPoints)
!
    ! Second derivatives
  do i = 1,nInter
    i0 = i*nPoints +1
    i1 = i0        +nPoints-1
!
    diffx_i(:,:) = -2.0D0*r(:,:,i)/l(i)
!
    do j = i,nInter
      j0 = j*nPoints+1
      j1 = j0+nPoints-1
!
      diffx_j(:,:)  =  2.0D0*r(:,:,j)/l(j)
!
    !   if differentiating twice on the same dimension
      if (i.eq.j) Then
       full_R(i0:i1,j0:j1) = matSder*diffx_j*diffx_i - matfder*(2.0D0/(l(i)*l(j)))
      else
       full_R(i0:i1,j0:j1) = matSder*diffx_j*diffx_i
      end if
    !   Writing the second derivatives in eq(2)
      if (i.ne.j) full_R(j0:j1,i0:i1) = transpose(Full_r(i0:i1,j0:j1))
    enddo
  enddo
!
!           Add constants to reflect the error in the energy and the
!           gradient, respectively.
!
  forall (j=1:nPoints) Full_R(j,j) = Full_R(j,j) + eps
  forall (j=nPoints+1:m_t) Full_R(j,j) = Full_R(j,j) + eps2
!
!           defining full_r has strictly positive define sec. 3 of
!           DOI: 10.1615/Int.J.UncertaintyQuantification.2013006809
  ! full_R = abs(full_R)
#ifdef _DEBUG_
  Call RecPrt('full_r Orig:','(14E10.2)',full_R,m_t,m_t)
#endif
!
  Call mma_deallocate(diffx_j)
  Call mma_deallocate(diffx_i)
  Call mma_deallocate(matFder)
  Call mma_deallocate(matSder)
  Call mma_deallocate(r)
  Call mma_deallocate(d)
END SUBROUTINE covarMatrix
