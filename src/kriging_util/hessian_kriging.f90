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
        use globvar
        Integer nInter,nPoints
        Real*8 x_(ndimx,1),ddy_(ndimx,ndimx)
!
!#define _Hess_Test
#ifdef _Hess_Test
        Real*8 Scale,Delta,Fact,tempX(ndimx,1),tgrad(ndimx),thgrad(ndimx)
        Real*8 HessT
!       Numerical Hessian threshold numHt from outside parameter
        HessT = numHt!1.0D-7
#endif
        nPoints = nPoints_save
        nInter = nInter_save
!
        ! write (6,*) 'nPoints', nPoints
        ! write (6,*) 'nInter', nInter
        npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
        nx(:,:) = x_
!
        ! Write(6,*) 'Entro hess'
! Analitical Hessian of GEK
        ! Call RecPrt('full_r:',  ' ',full_R,m_t,m_t)
        ! write(6,*) 'Kv:',Kv
        ! if (anHe) then
        call covarvector(2,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        call predict(2,nPoints,nInter)
        ddy_ = hpred(npx,:,:)
!               Factor to compare with Hessian Roland
                ! do i = 1,nInter
                !         hpred(npx,i,i) = 2.0D0*hpred(npx,i,i)
                ! enddo
        ! else
#ifdef _Hess_Test
! Numerical Hessian of GEK
        write(6,*) 'Begining Numerical Hessian'
!
        hpred = 0
        Scale=0.01D0
        write(6,*) 'Hess Threshold',numHt
        ! Call RecPrt('x = qInt',' ',x,nInter,nPoints)
        do i = 1,nInter
                tempx = x_
                ! Call RecPrt('x = qInt',' ',x_,nInter,1)
                Delta = 1.0D-5!Max(Abs(x_(i,1)),1.0D-5)*Scale
                ! write (6,*) 'x-qsave', i, x_(i,1)
                ! write (6,*) 'Delta', Delta
                ! write (6,*) 'nx', nx
                tempx(i,1) = x_(i,1) + Delta
                ! Call RecPrt('x(1,i) + Delta = qIntp',' ',tempx(:,1),nInter,1)
!
                Call Gradient_Kriging(tempx(:,1),tgrad,ndimx)
!
                tempx(i,1) = x_(i,1) - Delta
                Call Gradient_Kriging(tempx(:,1),thgrad,ndimx)
!
                ! Call RecPrt('x(1,i) - Delta = qIntm',' ',tempx(:,1),nInter,1)
                ! Call RecPrt('tgrad-dqp',' ',tgrad,nInter,1)
                ! Call RecPrt('thgrad-dqm',' ',thgrad,nInter,1)
                do j=1,nInter
                        Fact = 0.5D0
                        If (i.eq.j) Fact = 1.0D0
                        ! write(6,*) 'Hess i j',i,j
                        ! write(6,*) 'Fact, Delta',Fact, Delta
                        ! write(6,*) 'Hessian(i,j)',hpred(npx,i,j)
                        ! write(6,*) 'Hessian(j,i)',hpred(npx,j,i)
                        hpred(npx,i,j) = hpred(npx,i,j) + Fact*(tgrad(j)-thgrad(j))/(2.0D0*Delta)
                        hpred(npx,j,i) = hpred(npx,j,i) + Fact*(tgrad(j)-thgrad(j))/(2.0D0*Delta)
                        ! write(6,*) 'tgrad(j)',tgrad(j)
                        ! write(6,*) 'thgrad(j,i)',thgrad(j)
                        ! write(6,*) 'Hessian(i,j)',hpred(npx,i,j)
                        ! write(6,*) 'Hessian(j,i)',hpred(npx,j,i)

                enddo
        ! write (6,*) 'hpred(npx,i,:): ', i, hpred(npx,i,:)
        enddo
        do i = 1,nInter
                do j = 1,nInter
                        write(6,*) 'i,j',i,j,HessT
                        write(6,*) 'hpred, ddy_',hpred(npx,i,j),ddy_(i,j)
                        if (abs(ddy_(i,j)-hpred(npx,i,j)).gt.HessT) then
                                Write(6,*) 'Error in entry',i,',',j,'of the hessian matrix'
                                Call RecPrt('Anna Hess',' ',ddy_,nInter,nInter)
                                Call RecPrt('Num Hess',' ',hpred,nInter,nInter)
                                Write(6,*) 'abs(ddy_(i,j)+ HessT)',abs(ddy_(i,j)+ HessT)
                                Write(6,*) 'abs(ddy_(i,j)- HessT)',abs(ddy_(i,j)- HessT)
                                Write(6,*) 'abs(hpred(npx,i,j))',abs(hpred(npx,i,j))
                                Call Abend()
                        endif
                enddo
        enddo
        if (anHe.eqv..False.) ddy_ = hpred(npx,:,:)
#endif
        !--------temp---------------------
        ! write(6,*) 'x: ', x
        ! write(6,*) 'y: ', y
        ! write(6,*) 'nx: ', nx
        ! write(6,*) 'dy: ', dy
        ! write(6,*) 'l: ', l
        ! write(6,*) 'Kv: ',Kv
        ! write(6,*) 'Kriging Hessian, Analitical?', anHe
        ! ! Call RecPrt('Hess',' ',ddy_,nInter,nInter)
        ! write(6,*) '-------------------Ana Hess'
        ! call covarvector(2,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        ! call predict(2,nPoints,nInter)
        ! Call RecPrt('Anna Hess',' ',hpred(npx,:,:),nInter,nInter)
!---------------------------------
!
        return
      end
