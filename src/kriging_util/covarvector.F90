
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

subroutine covarVector(gh)

use kriging_mod, only: cv, cvMatFder, cvMatSder, cvMatTder, dl, Index_PGEK, l, nD, nInter, nInter_Eff, nPoints, rl, x, x0
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Three
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: gh
integer(kind=iwp) :: i, i0, i1, j, j0, j1, k, k0, k1, i_eff, j_eff, k_eff
real(kind=wp) :: sdiffxi, sdiffxj, sdiffxk
real(kind=wp), allocatable :: diffxi(:), diffxj(:), diffxk(:)

call mma_Allocate(diffxi,nPoints,label='diffxi')
call mma_Allocate(diffxj,nPoints,label='diffxj')
call mma_Allocate(diffxk,nPoints,label='diffxk')

cv = 0
i0 = 0
call defdlrl()

! Covariant Vector in kriging - First part of eq (4) in ref.
if (gh == 0) then

  call matern(dl,cv(1:nPoints,1,1),nPoints,1)
  call matderiv(1,dl,cvMatFDer,nPoints,1)
  do i_eff=1,nInter_eff
    i = Index_PGEK(i_eff)
    ! 1st derivatives second part of eq. (4)
    diffxi(:) = Two*rl(:,i)/l(i)
    i0 = nPoints+1+(i_eff-1)*(nPoints-nD)
    i1 = i0+(nPoints-nD)-1
    cv(i0:i1,1,1) = cvMatFder(1+nD:nPoints)*diffxi(1+nD:nPoints)
  end do
# ifdef _DEBUGPRINT_
  call RecPrt(' The covector for energies','(12(2x,E9.3))',cv(:,1,1),m_t,1)
# endif

! Covariant vector in Gradient Enhanced Kriging
else if (gh == 1) then

  call matderiv(1,dl,cvMatFder,nPoints,1)
  call matderiv(2,dl,cvMatSder,nPoints,1)
  do i=1,nInter
    diffxi(:) = Two*rl(:,i)/l(i)
    cv(1:nPoints,i,1) = -cvMatFder(1:nPoints)*diffxi(1:nPoints)
    do j_eff=1,nInter_eff
      j = Index_PGEK(j_eff)
      j0 = nPoints+1+(j_eff-1)*(nPoints-nD)
      j1 = j0+(nPoints-nD)-1
      diffxj(:) = -Two*rl(:,j)/l(j)
      if (i == j) then
        cv(j0:j1,i,1) = cvMatSder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)-cvMatFder(1+nD:nPoints)*(2/(l(i)*l(j)))
      else
        cv(j0:j1,i,1) = cvMatSder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)
      end if
    end do
  end do
# ifdef _DEBUGPRINT_
  call RecPrt(' The covector for gradients','(12(2x,E9.3))',cv(:,:,1),m_t,nInter)
# endif

else if (gh == 2) then

  ! print *,'covar vector calling deriv(3) for Kriging Hessian'
  call matderiv(1,dl,cvMatFder,nPoints,1)
  call matderiv(2,dl,cvMatSder,nPoints,1)
  call matderiv(3,dl,cvMatTder,nPoints,1)
  do i=1,nInter
    diffxi(:) = Two*rl(:,i)/l(i)
    sdiffxi = Two/l(i)**2
    do j=1,nInter
      diffxj(:) = Two*rl(:,j)/l(j)
      sdiffxj = Two/l(j)**2
      if (i == j) then
        cv(1:nPoints,i,j) = cvMatSder(:)*diffxi(:)*diffxj(:)+cvMatFder(:)*Two/(l(i)*l(j))
      else
        cv(1:nPoints,i,j) = cvMatSder(:)*diffxi(:)*diffxj(:)
      end if
      do k_eff=1,nInter_eff
        k = Index_PGEK(k_eff)
        diffxk(:) = Two*rl(:,k)/l(k)
        sdiffxk = Two/l(k)**2
        k0 = nPoints+1+(k_eff-1)*(nPoints-nD)
        k1 = k0+(nPoints-nD)-1
        if (i == j .and. j == k) then
          cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints)+ &
                          Three*cvMatSder(1+nD:nPoints)*diffxi(1+nD:nPoints)*sdiffxj
        else if (i == j) then
          cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints)+ &
                          cvMatSder(1+nD:nPoints)*diffxk(1+nD:nPoints)*sdiffxi
        else if (i == k) then
          cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints)+ &
                          cvMatSder(1+nD:nPoints)*diffxj(1+nD:nPoints)*sdiffxi
        else if (j == k) then
          cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints)+ &
                          cvMatSder(1+nD:nPoints)*diffxi(1+nD:nPoints)*sdiffxk
        else
          cv(k0:k1,i,j) = cvMatTder(1+nD:nPoints)*diffxi(1+nD:nPoints)*diffxj(1+nD:nPoints)*diffxk(1+nD:nPoints)
        end if
      end do
    end do
  end do
else
  write(u6,*) ' Illegal value of gh:',gh
  call Abend()
end if

call mma_deallocate(diffxi)
call mma_deallocate(diffxj)
call mma_deallocate(diffxk)

contains

subroutine defdlrl()

  use Constants, only: Zero
  use Definitions, only: iwp

  implicit none
  integer(kind=iwp) :: i, j

  dl(:) = Zero
  do i=1,nInter
    do j=1,nPoints
      rl(j,i) = (x(i,j)-x0(i))/l(i)
    end do
    dl(:) = dl(:)+rl(:,i)**2
  end do

end subroutine defdlrl

end subroutine covarvector
