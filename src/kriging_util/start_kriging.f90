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

      Subroutine Start_Kriging(iter,nInter,qInt,Grad,Energy,anAIT,pAIT,lbAI,npxAIT)
!#include "stdalloc.fh"
        use globvar
        Integer npxAIT,nInter,iter
        Real*8 qInt(nInter,iter+1),Grad(nInter,iter),Energy(iter),pAIT,lbAI(3)
        Logical anAIT
        ! if (allocated(x)) then
        !   deallocate (x,y,lb,dy,nx)
        ! else
          allocate (x(nInter,iter),y(iter),lb(3),dy(nInter*iter), &
                    nx(nInter,1))
        ! endif
        !iter=iterAI
        !nInter=nInterT
        npxAI=npxAIT
        pAI=pAIT
        anAI=anAIT
        npx=1 !npxAI
        nx=qInt(:,iter+1:iter+1)
        Write (6,*) 'Kriging values in Start Kriging'
        !Write (6,*) 'iter:', iter
        Write (6,*) 'nInter', nInter
        Write (6,*) 'npxAI', npxAI
        Write (6,*) 'Grad: ',Grad
        Write (6,*) 'Grad shape: ',shape(Grad)
        Write (6,*) 'Energy: ',Energy
        Write (6,*) 'npx: ',npx
        !---------------Testing old data

        !--------------------------------
        !anamat = anAI
        !p = pAI
        !NS = iter
        !dims = nInter
        !parameters comming from the predicted points
        !nx=x
        lb=lbAI
        !----for test puorposes
        lb(1)=0.00001
        lb(2)=500
        lb(3)=1000
        ! npx=4
        ! do i=1,int(npx)
        !   nx(1,i)=(real(i)-1.0)*4.0/real(npx-1)
        ! enddo
        !-----------
        x = qInt(1:nInter,1:iter)
        Write (6,*) 'nx: ',nx
        Write (6,*) 'nx size: ',size(nx)
        Write (6,*) 'nx shape: ',shape(nx)
        Write (6,*) 'x: ',x
        Write (6,*) 'x size: ',size(x)
        Write (6,*) 'x shape: ',shape(x)
        y = Energy
        Write (6,*) 'y: ',y
        Write (6,*) 'y size: ',size(y)
        Write (6,*) 'y shape: ',shape(y)
        do i=1,nInter
          do j=1,iter
            dy(j+(i-1)*iter) = Grad(i,j)
          enddo
        enddo
        Write (6,*) 'dy: ',dy
        Write (6,*) 'dy size: ',size(dy)
        Write (6,*) 'dy shape: ',shape(dy)
        call kernels(iter,nInter)
        Energy(iter+1)=pred(npx)
        Grad(:,iter+1)=gpred
        write(6,*) 'New values of Energy and grad', pred(npx), gpred
        deallocate (x,y,lb,dy,nx,l)
        !write (6,*) 'deallo x'
        deallocate (full_R,rl,dl,mat,Iden)
        !write (6,*) 'deallo full_R'
        deallocate (kv,pred,gpred,hpred,var,sigma,cv,ll)
        !write (6,*) 'deallo kv,pred,gpred,hpred,var,sigma,cv'
        return
      end
