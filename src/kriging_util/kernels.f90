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
            integer i,z,j,iter,nInter,lm!,lm
            ! real*8 tpgh(3,int(lb(3)))
!
            call miden(iter)
            z=int(lb(3))
!
            ! If (make_parameters) Then
!
                ! do j = 1,nInter
                !     do i = 1,z
                !      l(j)=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                !     enddo
                ! enddo
!To be change for the optmization of the l's (the right width of the Mat'ern function)
            do i = 1,z
!In this particullary case the l(j) it does not depend on the dimensionality
!It is the same.
!Could be the case that depends on the dimesionality then for every dimension the
!l changes
                do j = 1,nInter
                    l(j)=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                    ! write (6,*) 'l:',l(j)
                enddo
                call covarmatrix(iter,nInter)
                !Call RecPrt('Full_R:',' ',full_r,m_t,m_t)
                call k(iter)
                ! write (6,*) 'lh, variance:',lh, variance
                ! call covarvector(0,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
                ! call predict(0,iter,nInter)
                ll(i)=lh
                ! tpgh(1,i)=pred(npx)
                ! call covarvector(1,iter,nInter)
                ! call predict(1,iter,nInter)
                ! ll(2,i)=lh
                ! tpgh(2,i)=pred(npx)
            enddo
!
            ! lm_save = MaxLoc(abs(ll(:)),dim=nInter)
            lm = MaxLoc(ll,dim=nInter)
            do j = 1,nInter
                l(j)=lb(1)+(lm-1)*(lb(2)-lb(1))/(lb(3)-1)
            enddo
            Call covarmatrix(iter,nInter)
            Call RecPrt('Full_R:',' ',full_r,m_t,m_t)
            call k(iter)
            ! pred(npx)=tpgh(1,lm)
            Write(6,*) 'Likelihood function pred index lm,l,ll(lm),variance',lm,ll(lm),l,variance
!           lm=MaxLoc(abs(ll(2,:)),dim=nInter)
            ! gpred(npx)=tpgh(2,lm)
            ! Write(6,*) 'Likelihood function gpred',lm,gpred(npx)
            ! make_parameters=.False.
            ! lm_save = lm
!
!             Else
! !
!                 i = lm_save
!                 do j = 1,nInter
!                     l(j)=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
!                 enddo
!                     call covarvector(0,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
!                     call predict(0,iter,nInter)
!                     ll(1,i)=lh
!                     tpgh(1,i)=pred(npx)
!                     call covarvector(1,iter,nInter)
!                     call predict(1,iter,nInter)
!                     ll(2,i)=lh
!                     tpgh(2,i)=pred(npx)
!             pred(npx)=tpgh(1,lm_save)
!             Write(6,*) 'Likelihood function pred',lm_save,pred(npx)
!             gpred(npx)=tpgh(2,lm_save)
!             Write(6,*) 'Likelihood function gpred',lm_save,gpred(npx)
!             make_parameters=.False.
! !
!             End If
        END

        subroutine miden(iter)
            use globvar
            integer j,iter
            iden=0
            forall(j=1:iter) iden(j,j)=1
        end subroutine
