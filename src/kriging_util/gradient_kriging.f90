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
Subroutine Gradient_Kriging(x_,dy_,ndimx)
  use kriging_mod
  Implicit None
  Integer nInter,nPoints, ndimx
  Real*8 x_(ndimx,1),dy_(ndimx)
!
!#define _Grad_Test
#ifdef _Grad_Test
  Integer i
  Real*8 Delta,tpred,thpred
  Real*8 GradT
  GradT = 1.0D-6
#endif
  nPoints=nPoints_save
  nInter=nInter_save
!
  npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
  nx(:,:) = x_
!
  ! Write(6,*) 'Entro grad'
  call covarvector(1) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(1,nPoints,nInter)
  dy_=gpred(npx,:)
!
#ifdef _Grad_Test
 ! Numerical Gradient of GEK
  write(6,*) 'Begining Numerical Gradient'
  Delta = 1.0D-4!Max(Abs(x_(i,1)),1.0D-5)*Scale
  do i = 1,nInter
    nx(:,:) = x_
!
    nx(i,1) = x_(i,1) + Delta
!
    call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
    call predict(0,nPoints,nInter)
    tpred = pred(npx)
!
    nx(i,1) = x_(i,1) - Delta
!
    call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
    call predict(0,nPoints,nInter)
    thpred = pred(npx)
!
    gpred(npx,i) = (tpred-thpred)/(2.0D0*Delta)
!
  enddo
  write(6,*) 'Gradient Threshold',GradT
  do i = 1,nInter
    write(6,*) 'i',i
    write(6,*) 'gpred, dy_',gpred(npx,i),dy_(i)
    if (abs(dy_(i)-gpred(npx,i)).gt.GradT) then
      Write(6,*) 'Error in entry',i,'of the gradient vector'
      Call RecPrt('Anna Grad',' ',dy_,1,nInter)
      Call RecPrt('Num Grad',' ',gpred,1,nInter)
      Write(6,*) 'abs(dy_(i,j)+ HessT)',abs(dy_(i)+ GradT)
      Write(6,*) 'abs(dy_(i,j)- HessT)',abs(dy_(i)- GradT)
      Write(6,*) 'abs(gpred(npx,i))',abs(gpred(npx,i))
      Call Abend()
    endif
  enddo
#endif
  return
End Subroutine Gradient_Kriging
