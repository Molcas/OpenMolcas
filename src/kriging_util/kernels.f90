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
            integer i,z
            real*8 tpred(npx)
            m_t=NS*(1+dims)
            allocate (full_R(m_t,m_t))
            allocate (y(NS),dy(NS),rl(NS,NS),dl(NS,NS),mat(NS,NS),Iden(NS,NS),nx(npx))
            allocate (kv(m_t),pred(npx),var(npx),sigma(npx),cv(m_t,npx))
            call miden()
            z=int(lb(3))
            do i = 1,z
                l=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                call covarmatrix()
                call k()
                call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                call predict()
            enddo
            tpred=pred
            call covarvector(1)
            call predict()
            gpred=pred
            call covarvector(2)
            call predict()
            hpred=pred
            pred=tpred
        END

        subroutine miden()
            use globvar
            integer j
            iden=0
            forall(j=1:ns) iden(j,j)=1
        end subroutine
