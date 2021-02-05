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

use kriging_mod

implicit none
integer ndimx
real*8 x0_(ndimx), dy_(ndimx)
!#define _Grad_Test
#ifdef _Grad_Test
integer i
real*8 Delta, tpred, thpred
real*8 GradT
GradT = 1.0d-6
#endif

!x0 is the n-dimensional vector of the coordinates at which the gradient is evaluated.
! subroutine
x0(:) = x0_(:)

! Write(6,*) 'Entro grad'
call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
call predict(1)
dy_(:) = gpred(:)

#ifdef _Grad_Test
! Numerical Gradient of GEK
write(6,*) 'Begining Numerical Gradient'
Delta = 1.0d-4!Max(Abs(x_(i,1)),1.0D-5)*Scale
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

  gpred(i) = (tpred-thpred)/(2.0d0*Delta)

end do
write(6,*) 'Gradient Threshold',GradT
do i=1,nInter
  write(6,*) 'i',i
  write(6,*) 'gpred, dy_',gpred(i),dy_(i)
  if (abs(dy_(i)-gpred(i)) > GradT) then
    write(6,*) 'Error in entry',i,'of the gradient vector'
    call RecPrt('Anna Grad',' ',dy_,1,nInter)
    call RecPrt('Num Grad',' ',gpred,1,nInter)
    write(6,*) 'abs(dy_(i,j)+ HessT)',abs(dy_(i)+GradT)
    write(6,*) 'abs(dy_(i,j)- HessT)',abs(dy_(i)-GradT)
    write(6,*) 'abs(gpred(i))',abs(gpred(i))
    call Abend()
  end if
end do
#endif

return

end subroutine Gradient_Kriging
