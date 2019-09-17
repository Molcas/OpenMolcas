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
        use globvar
        Integer nInter,nPoints
        Real*8 x_(ndimx,1),dy_(ndimx)
!
!#define _Grad_Test
#ifdef _Grad_Test
        Real*8 Scale,Delta,Fact,tempX(ndimx,1),tgrad(ndimx),thgrad(ndimx)
        Real*8 GradT
        !       Numerical Gradient-Hessian threshold numHt from outside parameter
        GradT = numHt!1.0D-7
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
        call covarvector(1,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        call predict(1,nPoints,nInter)
        dy_=gpred(npx,:)
!
#ifdef _Grad_Test
       ! Numerical Gradient of GEK
        do i = 1,nInter
!
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
        return
      end
