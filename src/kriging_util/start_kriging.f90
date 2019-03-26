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

      Subroutine Start_Kriging(iterT,nInterT,qInt,Grad,Energy,anAIT,pAIT,lbAI,npxAIT)
        use globvar
        Integer npxAIT,nInterT,iterT
        Real*8 qInt(nInterT,iterT+1),Grad(nInterT,iterT),Energy(iterT),pAIT,lbAI(3)
        Logical anAIT
        allocate (x(nInterT,iterT),y(iterT),lb(3),dy(nInterT*iterT),nx(nInterT,1),l(nInterT))
        iter=iterT
        nInter=nInterT
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
        Write (6,*) 'Grad size: ',size(Grad)
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
        ! do i=1,int(npx)
        !   nx(1,i)=(real(i)-1.0)*4.0/real(npx-1)
        ! enddo
        lb=lbAI
        lb(1)=0.1
        lb(2)=200
        lb(3)=500
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
        call kernels()
      end