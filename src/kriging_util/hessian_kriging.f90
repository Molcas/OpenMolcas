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
        Real*8 x_(ndimx,1),ddy_(ndimx,ndimx),tgrad(ndimx),thgrad(ndimx)
        Real*8 Scale,Delta,Fact,tempX(nInter_save,nPoints_save)
!
        nPoints=nPoints_save
        nInter=nInter_save
!
        write (6,*) 'nPoints', nPoints
        write (6,*) 'nInter', nInter
        npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
        nx(:,:) = x_
!
        ! Write(6,*) 'Entro hess'
! Analitical Hessian of GEK
        ! Call RecPrt('full_r:',  ' ',full_R,m_t,m_t)
        ! write(6,*) 'Kv:',Kv
        if (anHe) then
                call covarvector(2,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                call predict(2,nPoints,nInter)
!               Factor to compare with Hessian Roland
                ! do i = 1,nInter
                !         hpred(npx,i,i) = 2.0D0*hpred(npx,i,i)
                ! enddo
        else
! Numerical Hessian of GEK
                ! write(6,*) 'Begining Numerical Hessian'
                ! wrire (6)
                hpred = 0
                Scale=0.01D0
                ! Call RecPrt('x = qInt',' ',x,nInter,nPoints)
                do i = 1,nInter
                        !qInt_Save= x(i,nPoints)
                        tempx = x
                        ! Call RecPrt('x = qInt',' ',x,nInter,nPoints)
                        Delta = Max(Abs(x(i,nPoints)),1.0D-5)*Scale
                        ! write (6,*) 'x-qsave', i, x(i,nPoints)
                        ! write (6,*) 'Delta', Delta
                        ! write (6,*) 'nx', nx
                        tempx(i,nPoints) = x(i,nPoints) + Delta
                        ! Call RecPrt('x(1,i) + Delta = qIntp',' ',tempx(:,nPoints),nInter,1)
                        Call Gradient_Kriging(tempx(:,nPoints),tgrad,ndimx)
                        !
                        tempx(i,nPoints) = x(i,nPoints) - Delta
                        ! Call RecPrt('x(1,i) - Delta = qIntm',' ',tempx(:,nPoints),nInter,1)
                        call Gradient_Kriging(tempx(:,nPoints),thgrad,ndimx)
                        ! Call RecPrt('tgrad-dqp',' ',tgrad,nInter,1)
                        ! Call RecPrt('thgrad-dqm',' ',thgrad,nInter,1)
                        do j=i,nInter
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
        endif
        ddy_=hpred(npx,:,:)
        !--------temp---------------------
        ! write(6,*) 'x: ', x
        ! write(6,*) 'y: ', y
        ! write(6,*) 'nx: ', nx
        ! write(6,*) 'dy: ', dy
        ! write(6,*) 'l: ', l
        ! write(6,*) 'Kv: ',Kv
        ! write(6,*) 'Kriging Hessian, Analitical?', anHe
        ! Call RecPrt('Hess',' ',ddy_,nInter,nInter)
        ! write(6,*) '-------------------Ana Hess'
        ! call covarvector(2,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        ! call predict(2,nPoints,nInter)
        ! Call RecPrt('Anna Hess',' ',hpred(npx,:,:),nInter,nInter)
        ! do i = 1,nInter
        !         hpred(npx,i,i) = 2.0*hpred(npx,i,i)
        ! enddo
        ! Call RecPrt('Anna Hess w factor',' ',hpred(npx,:,:),nInter,nInter)
!---------------------------------
!
        return
      end
