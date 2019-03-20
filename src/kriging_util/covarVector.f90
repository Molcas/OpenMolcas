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
        SUBROUTINE covarVector(gh)
            use globvar
            integer i,i0,i1,j,gh,k
            real*8 tmat(iter,npx),tmat2(iter,npx)!,c(m_t,npx)
            deallocate (rl,dl,mat)
            allocate (rl(iter,npx),dl(iter,npx),mat(iter,npx))
            cv=0
            tmat=0
            tmat2=0
            dl=0
            ! Covariant Vector in kriging
            if (gh.eq.0) then
                do i=1,nInter
                    do j=1,iter
                        do k=1,int(npx)
                            rl(j,k)=(x(i,j)-nx(k))/l(i)
                        enddo
                    enddo
                    !write(6,*) 'CV-rl',i,rl
                    dl=dl+rl**2
                    !write(6,*) 'CV-dl',i,dl
                    !write(6,*) 'CV-Mat',i,mat
                    !tmat=tmat+mat
                enddo
                call matern(dl,size(dl,1),size(dl,2))
                cv(1:iter,:,1)=mat
            endif
            ! Covariant vector in Gradient Enhanced Kriging
            do i=1,nInter
                do j=1,iter
                    do k=1,int(npx)
                        rl(j,k)=(x(i,j)-nx(k))/l(i)
                    enddo
                enddo
                dl=rl**2
                ! if (gh.eq.0) then
                !     call matern(dl,size(dl,1),size(dl,2))
                !     cv(1:iter,:,i)=mat
                ! endif
                if (gh.eq.0) then
                    i0=i*iter+1
                    i1=i0+iter-1
                else
                    i0=(i-1)*iter+1
                    i1=i0+iter-1
                endif
!                print *,'covar vector calling deriv(1)'
                call matderiv(1,size(dl,1),size(dl,2))
                tmat=mat
                mat=mat*(2*rl)/l(i)
                !Write (6,*) 'CV - Deriv: ',mat
                if(gh.ge.1) then
!                    print *,'covar vector calling deriv(2) for Kriging Gradients'
                    if (gh.eq.1) then
                        cv(i0:i1,:,i)=-mat
                        i0=i*iter+1
                        i1=i0+iter-1
                    endif
                    call matderiv(2,size(dl,1),size(dl,2))
                    tmat2=mat
                    mat=mat*(-4*rl**2/l(i)**2)+tmat*(-2/l(i)**2)
                    Write (6,*) 'CV - Krig Grad: ',mat
                endif
                if(gh.eq.2) then
!                    print *,'covar vector calling deriv(3)for Kriging Hessian'
                    cv(i0:i1,:,i)=mat
                    i0=i*iter+1
                    i1=i0+iter-1
                    call matderiv(3,size(dl,1),size(dl,2))
                    mat=mat*8*rl**4/l(i)**3+tmat2*(16*rl**2/l(i)**3 + 4/l(i)**2)+tmat*(4/l(i)**3)
                    Write (6,*) 'CV - Krig Hessian: ',mat
                endif
                cv(i0:i1,:,i)=mat
            enddo
            ! write (6,*) 'CV: ',CV
            ! Write (6,*) 'CV size: ',size(CV)
            ! Write (6,*) 'CV shape: ',shape(CV)
        END
