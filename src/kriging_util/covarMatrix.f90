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
        SUBROUTINE covarMatrix(iter,nInter)
            use globvar
            integer i,j,i0,i1,j0,j1,k,kl,iter,nInter
            real*8 tmat(iter,iter),diffx(iter,iter),diffx0(iter,iter)!,c(iter,iter)
            deallocate (rl,dl,mat)
            allocate (rl(iter,iter),dl(iter,iter), &
                     mat(iter,iter))
            ! Write (6,*) 'Covar Matrix'
            full_R=0
            !c=exp(-sqrt((2.0*p+1.0)*dl))
            tmat=0
            dl=0
            ! Covariant Matrix in kriging
            do i=1,nInter
                do k=1,iter
                    do kl=1,iter
                        rl(k,kl)=(x(i,k)-x(i,kl))/l(i)
                    end do
                end do
                dl=dl+rl**2
            end do
            call matern(dl,iter,iter)
            full_R(1:iter,1:iter)=mat+iden*eps
            ! Covariant matrix in Gradient Enhanced Kriging
            do i=1,nInter
                do k=1,iter
                    do kl=1,iter
                        rl(k,kl)=(x(i,k)-x(i,kl))/l(i)
                    end do
                end do
                dl=rl**2
                diffx=2.0*rl/l(i)
                diffx0=-2.0*rl/l(i)
                i0=i*iter+1
                i1=i0+iter-1
                do j=i,nInter+1
                    j0=(j-1)*iter+1
                    j1=j0+iter-1
                    if(j.eq.1) then
                        call matderiv(1,iter,iter)
                        tmat=mat
                        mat=mat*diffx
                    else
                        call matderiv(2,iter,iter)
                        mat=mat*diffx*diffx0-tmat*(2/l(i)**2)
                    endif
                    full_R(i0:i1,j0:j1)=mat
                    if (j.ne.nInter+1) then
                        full_R(j0:j1,i0:i1)=transpose(mat)
                    else
                        full_R(i0:i1,j0:j1)=full_R(i0:i1,j0:j1)+iden*eps
                    endif
                enddo
            enddo
            ! write (6,*) 'full_R: ',full_R
            ! Write (6,*) 'full_R size: ',size(full_R)
            ! Write (6,*) 'full_R shape: ',shape(full_R)
        END