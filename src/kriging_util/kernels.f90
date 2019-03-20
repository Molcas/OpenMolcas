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
        SUBROUTINE kernels()
            use globvar
            integer i,z,j
            real*8 tpred(npx)
            m_t=iter*(1+nInter)
            allocate (full_R(m_t,m_t))
            allocate (rl(iter,iter),dl(iter,iter), &
                      mat(iter,iter),Iden(iter,iter))
            allocate (kv(m_t),pred(npx),var(npx),sigma(npx),cv(m_t,npx,nInter))
            call miden()
            z=int(lb(3))
            Write (6,*) 'Kernels l', z
            do j = 1,nInter
                do i = 1,z
                    l((j))=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                    call covarmatrix()
                    call k()
                    call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                    call predict(0)
                    call covarvector(1)
                    call predict(1)
                    call covarvector(2)
                    call predict(2)
                enddo
            enddo
            tpred=pred
            ! call covarvector(1)
            ! call predict(1)
            gpred=pred
            ! call covarvector(2)
            ! call predict(1)
            hpred=pred
            pred=tpred
        END

        subroutine miden()
            use globvar
            integer j
            iden=0
            forall(j=1:iter) iden(j,j)=1
        end subroutine