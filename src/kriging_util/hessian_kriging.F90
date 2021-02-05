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

use kriging_mod

implicit none
integer ndimx
real*8 x0_(ndimx), ddy_(ndimx,ndimx)
!#define _Hess_Test
#ifdef _Hess_Test
real*8 Scale, Delta, Fact, tgrad(ndimx), thgrad(ndimx)
real*8 HessT, tmp
integer i, j
HessT = 1.0d-3
#endif

!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
x0(:) = x0_(:)

call covarvector(2) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
call predict(2)
ddy_(:,:) = hpred(:,:)

#ifdef _Hess_Test
! Numerical Hessian of GEK
write(6,*) 'Begining Numerical Hessian'

hpred(:,:) = 0.0d0
Scale = 0.01d0
write(6,*) 'Hess Threshold',HessT

do i=1,nInter
  tmp = x0(i)

  Delta = 1.0d-5!Max(Abs(x_(i,1)),1.0D-5)*Scale

  x0(i) = tmp+Delta
  call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(1)
  tgrad = gpred(:)

  x0(i) = tmp-Delta
  call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(1)
  thgrad = gpred(:)

  do j=1,nInter
    hpred(i,j) = (tgrad(j)-thgrad(j))/(2.0d0*Delta)
  end do
  x0(i) = tmp
end do
! Comparing Analytical solution with Numerical
do i=1,nInter
  do j=1,nInter
    write(6,*) 'i,j',i,j
    write(6,*) 'hpred, ddy_',hpred(i,j),ddy_(i,j)
    if (abs(ddy_(i,j)-hpred(i,j)) > HessT) then
      write(6,*) 'Error in entry',i,',',j,'of the hessian matrix'
      call RecPrt('Anal Hess',' ',ddy_,nInter,nInter)
      call RecPrt('Num Hess',' ',hpred,nInter,nInter)
      write(6,*) 'abs(ddy_(i,j)+ HessT)',abs(ddy_(i,j)+HessT)
      write(6,*) 'abs(ddy_(i,j)- HessT)',abs(ddy_(i,j)-HessT)
      write(6,*) 'abs(hpred(i,j))',abs(hpred(i,j))
      call Abend()
    end if
  end do
end do
#endif

return

end subroutine Hessian_Kriging
