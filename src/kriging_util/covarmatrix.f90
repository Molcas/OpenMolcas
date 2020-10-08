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
SUBROUTINE covarMatrix()
  use kriging_mod
  Implicit None
#include "stdalloc.fh"
  integer i,j,i0,i1,j0,j1,k,nInter
  Real*8, Allocatable :: diffx_j(:,:), diffx_i(:,:), matFder(:,:),&
                         matSder(:,:), r(:,:,:), d(:,:)
!#define _DEBUG_
nInter=nInter_save

!**********************************************************************
!
! Allocate temporary memory
!
  Call mma_Allocate(diffx_j,nPoints,nPoints,Label="diffx_j")
  Call mma_Allocate(diffx_i,nPoints,nPoints,Label="diffx_i")
  Call mma_Allocate(matFder,nPoints,nPoints,Label="matFder")
  Call mma_Allocate(matSder,nPoints,nPoints,Label="matSder")
  Call mma_Allocate(r,nPoints,nPoints,nInter,Label="r")
  Call mma_Allocate(d,nPoints,nPoints,Label="d")
!
!**********************************************************************
!
! Compute the distances between the sample points using
! the characteristic length.
!
  full_R(:,:) = 0.0D0
  d(:,:) = 0.0D0
  diffx_j(:,:) = 0.0D0
  diffx_i(:,:) = 0
!
  do i=1,nInter

    do k=1,nPoints
      do j=1,nPoints
        r(k,j,i)=(x(i,k)-x(i,j))/l(i)
      end do
    end do

!   Accumulate contributions to the square of the individual distances.

    d(:,:) = d(:,:) + r(:,:,i)**2

#ifdef _DEBUG_
    Call RecPrt('r',' ',r(1,1,i),nPoints,nPoints)
#endif

  end do
!
!**********************************************************************
!
#ifdef _DEBUG_
  Call RecPrt('l',' ',l,1,nInter)
  Call RecPrt('x',' ',x,nInter,nPoints)
  Call RecPrt('d',' ',d,nPoints,nPoints)
#endif
!
!**********************************************************************
!**********************************************************************
!
! Now evaluate the covariance function over all the distances. For GEK
! we will need gradients and 2nd order derivatives of the covariance
! function too.
!
! Currently we use the 5/2 Matern function as the covariance function
!
!**********************************************************************
!
! Note that we will evaluate the derivative of the covariance function
! w.r.t d. For the full derivative this has to be complemented by
! the derivative of d w.r.t the individual components of the coordinates.
!
!**********************************************************************
! 1) Evaluate the covariance function for all the distances.
!
  Call matern(d, full_R(1:nPoints,1:nPoints), nPoints, nPoints)
!
! Writing the covariant matrix in GEK (eq 2 of doi:10.1007/s00366-015-0397)
!
!**********************************************************************
!
! 2) Evaluate first derivatives of the covariance function with respect to d at all distances.
!
  Call matderiv(1, d, MatFder, nPoints, nPoints)
!
! Covariant matrix in Gradient Enhanced Kriging (eq 2 of doi:10.1007/s00366-015-0397):
!
! First line and first column derivative in Psi matrix
!
  do i=1,nInter      ! Loop over component of the coordinate to differentiate
!
!   Compute the range of the block in the covariance matrix.
!
    i0 = nPoints + 1 + (i-1)*(nPoints-nD)
    i1 = i0 + (nPoints-nD) - 1
!
!   Do an on-the-fly evaluation of the dervative w.r.t x_i
    diffx_i(1:nPoints,1+nD:nPoints) = -2.0D0*r(1:nPoints,1+nD:nPoints,i)/l(i)
!
!   Writing the 1st row of 1st derivatives with respect the coordinates
!
    full_R(1:nPoints,i0:i1) = matFDer(1:nPoints,1+nD:nPoints)  &
                            * diffx_i(1:nPoints,1+nD:nPoints)

  enddo
! Complete by filling in the transpose blocks

  full_R(nPoints+1:m_t,1:nPoints) = Transpose(Full_R(1:nPoints,nPoints+1:m_t))
!
!**********************************************************************
!
! 3) Evaluate the second derivatives.
!
! Matern second derivative with respect to d
!
  call matderiv(2, d, matSder, nPoints, nPoints)
!
    ! Second derivatives
  do i = 1,nInter
    i0 = nPoints + 1 + (i-1)*(nPoints-nD)
    i1 = i0 + (nPoints-nD) - 1
!
    diffx_i(1+nD:nPoints,1+nD:nPoints) = -2.0D0*r(1+nD:nPoints,1+nD:nPoints,i)/l(i)
!
    do j = i,nInter
      j0 = nPoints + 1 + (j-1)*(nPoints-nD)
      j1 = j0 + (nPoints-nD) - 1
!
      diffx_j(1+nD:nPoints,1+nD:nPoints)  =  2.0D0*r(1+nD:nPoints,1+nD:nPoints,j)/l(j)
!
    !   if differentiating twice on the same dimension
      full_R(i0:i1,j0:j1) = matSder(1+nD:nPoints,1+nD:nPoints) &
                          * diffx_j(1+nD:nPoints,1+nD:nPoints) &
                          * diffx_i(1+nD:nPoints,1+nD:nPoints)

      if (i.eq.j) full_R(i0:i1,j0:j1) = full_R(i0:i1,j0:j1) &
                                      - matFder(1+nD:nPoints,1+nD:nPoints)*(2.0D0/(l(i)*l(j)))

    !   Writing the second derivatives in eq(2)
      if (i.ne.j) full_R(j0:j1,i0:i1) = transpose(Full_r(i0:i1,j0:j1))
    enddo

  enddo
!
!           Add constants to reflect the error in the energy and the
!           gradient, respectively.
!
  do j=1,m_t
    if (j.le.nPoints) then
      Full_R(j,j) = Full_R(j,j) + eps
    else
      Full_R(j,j) = Full_R(j,j) + eps2
    end if
  end do
!
!           defining full_r has strictly positive define sec. 3 of
!           doi:10.1615/Int.J.UncertaintyQuantification.2013006809
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
