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
        Real*8 qInt(nInterT,iterT),Grad(nInterT,iterT),Energy(iterT),pAIT,lbAI(3)
        Logical anAIT
        allocate (x(nInterT,iterT),y(iterT),lb(3),dy(nInterT*iterT),nx(iterT),l(nInterT))
        iter=iterT
        nInter=nInterT
        npxAI=npxAIT
        pAI=pAIT
        anAI=anAIT
        Write (6,*) 'Kriging values in Start Kriging'
        Write (6,*) 'iter:', iter
        Write (6,*) 'nInter', nInter
        Write (6,*) 'npxAI', npxAI
        Write (6,*) 'Grad: ',Grad
        Write (6,*) 'Grad size: ',size(Grad)
        Write (6,*) 'Grad shape: ',shape(Grad)
        Write (6,*) 'Energy: ',Energy
        !---------------Testing old data

        !--------------------------------
        !anamat = anAI
        !p = pAI
        !NS = iter
        !dims = nInter
        npx=iter !npxAI
        do i=1,int(npx)
          nx(i)=(real(i)-1.0)*4.0/real(npx-1)
        enddo
        lb=lbAI
        x = qInt
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