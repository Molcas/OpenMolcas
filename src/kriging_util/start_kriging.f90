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

      Subroutine Start_Kriging(iter,nInter,qInt,Grad,Energy,pAIT,lbAI)
        use globvar
        Integer nInter,iter
        Real*8 qInt(nInter,iter+1),Grad(nInter,iter),Energy(iter),pAIT,lbAI(3)
!
        allocate (x(nInter,iter),y(iter),lb(3),dy(nInter*iter), &
                    nx(nInter,1))
!
        pAI=pAIT
        npx=1 !npxAI
        nx=qInt(:,iter+1:iter+1)
        lb=lbAI
        x = qInt(1:nInter,1:iter)
        y = Energy
        do i=1,nInter
          do j=1,iter
            dy(j+(i-1)*iter) = Grad(i,j)
          enddo
        enddo
!
        call kernels(iter,nInter)
!
        Energy(iter+1)=pred(npx)
        Grad(:,iter+1)=gpred

        write(6,*) 'New values of Energy and grad', pred(npx), gpred
        deallocate (x,y,lb,dy,nx,l)
        deallocate (full_R,rl,dl,mat,Iden)
        deallocate (kv,pred,gpred,hpred,var,sigma,cv,ll)
!
        return
      end
