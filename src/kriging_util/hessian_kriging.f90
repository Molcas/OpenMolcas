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
Subroutine Hessian_Kriging(x_,ddy_,ndimx)
  use kriging_mod
  Implicit None
  Integer nInter,nPoints,ndimx
  Real*8 x_(ndimx,1),ddy_(ndimx,ndimx)
!
!#define _Hess_Test
#ifdef _Hess_Test
  Real*8 Scale,Delta,Fact,tgrad(ndimx),thgrad(ndimx)
  Real*8 HessT, tmp
  Integer i, j
  HessT = 1.0D-3
#endif
  nPoints = nPoints_save
  nInter = nInter_save
!
  npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
  nx(:,:) = x_
!
  call covarvector(2,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(2,nPoints,nInter)
  ddy_ = hpred(npx,:,:)
!
#ifdef _Hess_Test
! Numerical Hessian of GEK
  write(6,*) 'Begining Numerical Hessian'
!
  hpred = 0
  Scale=0.01D0
  write(6,*) 'Hess Threshold',HessT
!
  do i = 1,nInter
    tmp=nx(i,1)
!
    Delta = 1.0D-5!Max(Abs(x_(i,1)),1.0D-5)*Scale
!
    nx(i,1) = tmp + Delta
    call covarvector(1,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
    call predict(1,nPoints,nInter)
    tgrad=gpred(npx,:)
!
    nx(i,1) = tmp - Delta
    call covarvector(1,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
    call predict(1,nPoints,nInter)
    thgrad=gpred(npx,:)
!
    do j=1,nInter
      hpred(npx,i,j) = (tgrad(j)-thgrad(j))/(2.0D0*Delta)
    enddo
    nx(i,1) = tmp
  enddo
! Comparing Analytical solution with Numerical
  do i = 1,nInter
    do j = 1,nInter
      write(6,*) 'i,j',i,j
      write(6,*) 'hpred, ddy_',hpred(npx,i,j),ddy_(i,j)
      if (abs(ddy_(i,j)-hpred(npx,i,j)).gt.HessT) then
        Write(6,*) 'Error in entry',i,',',j,'of the hessian matrix'
        Call RecPrt('Anal Hess',' ',ddy_,nInter,nInter)
        Call RecPrt('Num Hess',' ',hpred,nInter,nInter)
        Write(6,*) 'abs(ddy_(i,j)+ HessT)',abs(ddy_(i,j)+ HessT)
        Write(6,*) 'abs(ddy_(i,j)- HessT)',abs(ddy_(i,j)- HessT)
        Write(6,*) 'abs(hpred(npx,i,j))',abs(hpred(npx,i,j))
        Call Abend()
      endif
    enddo
  enddo
#endif
!
  return
End Subroutine Hessian_Kriging
