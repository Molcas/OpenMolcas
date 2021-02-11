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

subroutine Gradient_Kriging(x0_,dy_,ndimx)

use kriging_mod, only: gpred, x0
use Definitions, only: wp, iwp

!#define _Grad_Test
#ifdef _Grad_Test
use Constants, only: Two, u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ndimx
real(kind=wp), intent(in) :: x0_(ndimx)
real(kind=wp), intent(out) :: dy_(ndimx)
#ifdef _Grad_Test
integer(kind=iwp) :: i
real(kind=wp) :: Delta, tpred, thpred
real(kind=wp), parameter :: GradT = 1.0e-6_wp
#endif

!x0 is the n-dimensional vector of the coordinates at which the gradient is evaluated.
! subroutine
x0(:) = x0_(:)

! Write(u6,*) 'Entro grad'
call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
call predict(1)
dy_(:) = gpred(:)

#ifdef _Grad_Test
! Numerical Gradient of GEK
write(u6,*) 'Begining Numerical Gradient'
Delta = 1.0e-4_wp !max(abs(x_(i,1)),1.0e-5_wp)*Scale
do i=1,nInter
  x0(:) = x0_(:)

  x0(i) = x0_(i)+Delta

  call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(0)
  tpred = pred

  x0(i) = x0_(i)-Delta

  call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(0)
  thpred = pred

  gpred(i) = (tpred-thpred)/(Two*Delta)

end do
write(u6,*) 'Gradient Threshold',GradT
do i=1,nInter
  write(u6,*) 'i',i
  write(u6,*) 'gpred, dy_',gpred(i),dy_(i)
  if (abs(dy_(i)-gpred(i)) > GradT) then
    write(u6,*) 'Error in entry',i,'of the gradient vector'
    call RecPrt('Anna Grad',' ',dy_,1,nInter)
    call RecPrt('Num Grad',' ',gpred,1,nInter)
    write(u6,*) 'abs(dy_(i,j)+ HessT)',abs(dy_(i)+GradT)
    write(u6,*) 'abs(dy_(i,j)- HessT)',abs(dy_(i)-GradT)
    write(u6,*) 'abs(gpred(i))',abs(gpred(i))
    call Abend()
  end if
end do
#endif

return

end subroutine Gradient_Kriging
