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
            real*8 tmat(NS,npx),tmat2(NS,npx)!,c(m_t,npx)
            deallocate (rl,dl,mat)
            allocate (rl(NS,npx),dl(NS,npx),mat(NS,npx))
            cv=0
            tmat=0
            tmat2=0
            ! Covariant Vector in kriging
            do i=1,dims
                do j=1,ns
                    do k=1,int(npx)
                        rl(j,k)=(x(i,j)-nx(k))/l
                    enddo
                enddo
                write(6,*) 'CV-rl',i,rl
                dl=rl**2
                write(6,*) 'CV-dl',i,dl
                call matern(dl,size(dl,1),size(dl,2))
                write(6,*) 'CV-Mat',i,mat
                tmat=tmat+mat
            enddo
            cv(1:NS,:)=tmat
            ! Covariant vector in Gradient Enhanced Kriging
            do i=1,dims
                do j=1,ns
                    do k=1,int(npx)
                        rl(j,k)=(x(i,j)-nx(k))/l
                    enddo
                enddo
                dl=rl**2
                i0=i*ns+1
                i1=i0+ns-1
!                print *,'covar vector calling deriv(1)'
                call matderiv(1,size(dl,1),size(dl,2))
                tmat=mat
                mat=mat*(2*rl)/l
                Write (6,*) 'CV - Deriv: ',mat
                if(gh.ge.1) then
!                    print *,'covar vector calling deriv(2) for Kriging Gradients'
                    call matderiv(2,size(dl,1),size(dl,2))
                    tmat2=mat
                    mat=mat*(-4*rl**2/l**2)+tmat*(-2/l**2)
                    Write (6,*) 'CV - Krig Grad: ',mat
                endif
                if(gh.eq.2) then
!                    print *,'covar vector calling deriv(3)for Kriging Hessian'
                    call matderiv(3,size(dl,1),size(dl,2))
                    mat=mat*8*rl**4/l**3+tmat2*(16*rl**2/l**3 + 4/l**2)+tmat*(4/l**3)
                    Write (6,*) 'CV - Krig Hessian: ',mat
                endif
                cv(i0:i1,:)=mat
            enddo
            write (6,*) 'CV: ',CV
            Write (6,*) 'CV size: ',size(CV)
            Write (6,*) 'CV shape: ',shape(CV)
        END
