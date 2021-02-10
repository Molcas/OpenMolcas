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

subroutine Hessian_Kriging(x0_,ddy_,ndimx)

use kriging_mod, only: hpred, x0
use Definitions, only: wp, iwp

!#define _Hess_Test
#ifdef _Hess_Test
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ndimx
real(kind=wp), intent(in) :: x0_(ndimx)
real(kind=wp), intent(out) :: ddy_(ndimx,ndimx)
#ifdef _Hess_Test
integer(kind=iwp) :: i, j
real(kind=wp) :: Delta, Fact, tmp
real(kind=wp), allocatable :: tgrad(:), thgrad(:)
real(kind=wp), parameter :: HessT = 1.0e-3_wp
#endif

!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
x0(:) = x0_(:)

call covarvector(2) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
call predict(2)
ddy_(:,:) = hpred(:,:)

#ifdef _Hess_Test
! Numerical Hessian of GEK
write(u6,*) 'Begining Numerical Hessian'

call mma_allocate(tgrad,ndimx,label='tgrad')
call mma_allocate(thgrad,ndimx,label='thgrad')

hpred(:,:) = Zero
write(u6,*) 'Hess Threshold',HessT

do i=1,nInter
  tmp = x0(i)

  Delta = 1.0e-5_wp !max(abs(x_(i,1)),1.0e-5_wp)*Scale

  x0(i) = tmp+Delta
  call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(1)
  tgrad = gpred(:)

  x0(i) = tmp-Delta
  call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(1)
  thgrad = gpred(:)

  do j=1,nInter
    hpred(i,j) = (tgrad(j)-thgrad(j))/(Two*Delta)
  end do
  x0(i) = tmp
end do
! Comparing Analytical solution with Numerical
do i=1,nInter
  do j=1,nInter
    write(u6,*) 'i,j',i,j
    write(u6,*) 'hpred, ddy_',hpred(i,j),ddy_(i,j)
    if (abs(ddy_(i,j)-hpred(i,j)) > HessT) then
      write(u6,*) 'Error in entry',i,',',j,'of the hessian matrix'
      call RecPrt('Anal Hess',' ',ddy_,nInter,nInter)
      call RecPrt('Num Hess',' ',hpred,nInter,nInter)
      write(u6,*) 'abs(ddy_(i,j)+ HessT)',abs(ddy_(i,j)+HessT)
      write(u6,*) 'abs(ddy_(i,j)- HessT)',abs(ddy_(i,j)-HessT)
      write(u6,*) 'abs(hpred(i,j))',abs(hpred(i,j))
      call Abend()
    end if
  end do
end do

call mma_deallocate(tgrad)
call mma_deallocate(thgrad)
#endif

return

end subroutine Hessian_Kriging
