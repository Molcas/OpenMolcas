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
      Subroutine Energy_Kriging(x_,y_,ndimx)
        use globvar
        Integer nInter,nPoints
        Real*8 x_(ndimx,1),y_
!
        nPoints=nPoints_save
        nInter=nInter_save
!
        npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
        nx = x_
!
!       Write(6,*) 'Entro predict x_', x_
!       Write(6,*) 'nx', nx
!       Write(6,*) 'ndimx', ndimx

        call covarvector(0,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        call predict(0,nPoints,nInter)
        y_=pred(npx)
!
        return
      end
!
!
      Subroutine Gradient_Kriging(x_,dy_,ndimx)
        use globvar
        Integer nInter,nPoints
        Real*8 x_(ndimx,1),dy_(ndimx)
!
        nPoints=nPoints_save
        nInter=nInter_save
!
        npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
        nx = x_
!
        ! Write(6,*) 'Entro grad'
        call covarvector(1,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        call predict(1,nPoints,nInter)
        dy_=gpred(npx,:)
!
        return
      end
!
      Subroutine Hessian_Kriging(x_,ddy_,ndimx)
        use globvar
        Integer nInter,nPoints
        Real*8 x_(ndimx,1),ddy_(ndimx,ndimx),tgrad(ndimx),thgrad(ndimx)
        ! Real*8 Scale,Delta
!
        nPoints=nPoints_save
        nInter=nInter_save
!
        npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
        nx = x_
!
        ! Write(6,*) 'Entro hess'
! Analitical Hessian of GEK
        if (anHe) then
                call covarvector(2,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                call predict(2,nPoints,nInter)
        else
! Numerical Hessian of GEK
                hpred = 0
                do i = 1,nInter
                        ! Scale=0.01D0
                        ! Delta = Max(Abs(nx(i,1)),1.0D-5)*Scale
                        call Gradient_Kriging(nx,tgrad,ndimx)
                        call Gradient_Kriging(nx+h,thgrad,ndimx)
                        do j=1,nInter
                                ! Fact = 0.5D0
                                ! If (i.eq.j) Fact = 1.0D0
                                hpred(npx,i,j) = hpred(npx,i,j) + (thgrad(j)-tgrad(j))/h
                                hpred(npx,j,i) = hpred(npx,i,j)
                        enddo
                        ! write (6,*) 'Delta: ', Delta
                enddo
        endif
        ddy_=hpred(npx,:,:)
        ! write(6,*) 'Kriging Hessian, Analitical?', anHe, ddy_
!--------temp
        call covarvector(2,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        call predict(2,nPoints,nInter)
        ! write(6,*) 'Kriging Hessian, Analitical', hpred(npx,:,:)
!
        return
      end
!
      Subroutine Dispersion_Kriging(x_,y_,ndimx)
        use globvar
        Integer nInter,nPoints
        Real*8 x_(ndimx,1),y_(npxAI)
!
        ! npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
!
        if (all(x_.ne.nx)) then
                nPoints=nPoints_save
                nInter=nInter_save
                nx = x_
                call covarvector(0,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                call predict(0,nPoints,nInter)
        endif
        y_ = sigma
        ! write (6,*) 'sigma', sigma
!
        return
      end
