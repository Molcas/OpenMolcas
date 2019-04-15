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
        SUBROUTINE kernels(iter,nInter)
            use globvar
            integer i,z,j,iter,nInter,lm
            real*8 tpgh(3,int(lb(3)))
!
            call miden(iter)
            z=int(lb(3))
!
            If (make_parameters) Then
!
                do j = 1,nInter
                    do i = 1,z
                     l(j)=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                    enddo
                enddo
                do i = 1,z
                    call covarmatrix(iter,nInter)
                    call k(iter)
                    call covarvector(0,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                    call predict(0,iter,nInter)
                    ll(1,i)=lh
                    tpgh(1,i)=pred(npx)
                    call covarvector(1,iter,nInter)
                    call predict(1,iter,nInter)
                    ll(2,i)=lh
                    tpgh(2,i)=pred(npx)
                enddo
                lm=MaxLoc(abs(ll(1,:)),dim=nInter)
                pred(npx)=tpgh(1,lm)
                Write(6,*) 'Likelihood function pred',lm,pred(npx)
    !           lm=MaxLoc(abs(ll(2,:)),dim=nInter)
                gpred(npx)=tpgh(2,lm)
                Write(6,*) 'Likelihood function gpred',lm,gpred(npx)
                make_parameters=.False.
                lm_save = lm
!
            Else
!
                i = lm_save
                do j = 1,nInter
                    l(j)=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                enddo
                    call covarvector(0,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                    call predict(0,iter,nInter)
                    ll(1,i)=lh
                    tpgh(1,i)=pred(npx)
                    call covarvector(1,iter,nInter)
                    call predict(1,iter,nInter)
                    ll(2,i)=lh
                    tpgh(2,i)=pred(npx)
            pred(npx)=tpgh(1,lm_save)
            Write(6,*) 'Likelihood function pred',lm_save,pred(npx)
            gpred(npx)=tpgh(2,lm_save)
            Write(6,*) 'Likelihood function gpred',lm_save,gpred(npx)
            make_parameters=.False.
!
            End If
        END

        subroutine miden(iter)
            use globvar
            integer j,iter
            iden=0
            forall(j=1:iter) iden(j,j)=1
        end subroutine
