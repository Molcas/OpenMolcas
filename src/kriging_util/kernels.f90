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
            ! real*8 temp_pred(npx,int(lb(3))),temp_gpred(npx,nInter,int(lb(3)))
!
            call miden(iter)
            z=int(lb(3))
!
            !temp nx for testing
        !   nx(1,1) = x(1,1)
        !   nx(2,1) = x(2,1)
        !   nx(3,1) = x(3,1)
        !   write (6,*) 'x',x
        !   write (6,*) 'y',y
        !   write (6,*) 'dy',dy
        !   write (6,*) 'nx',nx
!To be change for the optmization of the l's (the right width of the Mat'ern function)
            do i = 1,z
!In this particullary case the l(j) it does not depend on the dimensionality
!It is the same.
!Could be the case that depends on the dimesionality then for every dimension the
!l changes
                do j = 1,nInter
                    l(j)=lb(1)+(i-1)*(lb(2)-lb(1))/(lb(3)-1)
                enddo
!
                call covarmatrix(iter,nInter)
                call k(iter)
                ll(i)=lh
                !------testing
              write (6,*) 'di i,l,lh:',i,l(1),ll(i)
            !   call covarvector(0,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
            !   call predict(0,iter,nInter)
            ! !   temp_pred(:,i)=pred
            !   write(6,*) 'pred',pred
            !   call covarvector(1,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
            !   call predict(1,iter,nInter)
            ! !   temp_gpred(:,:,i)=gpred
            !   write(6,*) 'gpred',gpred
            !   call covarvector(2,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
            !   call predict(2,iter,nInter)
            !   !temp_hpred(:,:,i)=gpred
            !   !write(6,*) 'hpred',hpred
            !   Call RecPrt('hpred',  ' ',hpred(npx,:,:),nInter,nInter)
            !   write (6,*) '------------------------'
            !----------
            enddo
!
            lm = MinLoc(ll,dim=nInter)
            do j = 1,nInter
                l(j)=lb(1)+(lm-1)*(lb(2)-lb(1))/(lb(3)-1)
            enddo
            Call covarmatrix(iter,nInter)
            Call k(iter)
            write (6,*) 'optimazed l, lh:',l(1),ll(lm)
            if (blaAI) then
                write (6,*) ''
                write (6,*) 'Baseline (Trend Function) has been added with: ', blavAI
                write (6,*) 'to the value: ',sb
                write (6,*) '(the ordinary value is: ', sbO, ')'
            else
                if (mblAI) then
                    write (6,*) ''
                    write (6,*) 'Baseline (Trend Function) changed to the maximum value of the energy: ', sb
                    write (6,*) '(the ordinary value is: ', sbO, ')'
                else
                    if (blAI) then
        !    Trend Function eq (7) in ref.
                        write (6,*) ''
                        write (6,*) 'Baseline (Trend Function) changed to value: ', sb
                        write (6,*) '(the ordinary value is: ', sbO, ')'
                    endif
                endif
            endif
            !--------testing---------
            ! nx(1,1) = 0.0000000000006019
            ! nx(2,1) = 1.5477663075629295
            ! nx(3,1) = 2.7651245913505389
            ! call covarvector(2,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
            ! call predict(2,iter,nInter)
            ! write(6,*) 'hpred',hpred
            ! Call RecPrt('Update_: hpred',' ',hpred(1,:,:),nInter,nInter)
            !-------------------------
            ! Write(6,*) 'optimazed lh: ',lm,ll(lm),l(1)
            ! write(6,*) 'all ll:',ll
            ! write(6,*) 'all Ener:',temp_pred
            ! write(6,*) 'all Grad:',temp_gpred
        END

        SUBROUTINE setlkriging(lv,iter,nInter)
            use globvar
            integer iter,nInter,i
            real*8 lv
            do i = 1,nInter
                l(i)=lv
            enddo
            call covarmatrix(iter,nInter)
            call k(iter)
        END

        subroutine miden(iter)
            use globvar
            integer j,iter
            iden=0
            forall(j=1:iter) iden(j,j)=1
        end subroutine
