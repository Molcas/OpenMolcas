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
            real*8 tmat(m_t,npx),tmat2(m_t,npx)!,c(m_t,npx)
            deallocate (rl,dl,mat)
            allocate (rl(NS,npx),dl(NS,npx),mat(NS,npx))
            ! do i=1,dims
            !     do j=1,ns
            !         do k=1,int(npx)
            !             rl(i,j)=(x(i)-nx(k))/l
            !         enddo
            !     enddo
            ! enddo
            ! dl=rl**2
            ! c=exp(-sqrt((2.0*p+1.0)*dl))
            Write (6,*) 'Covar Vector'
            cv=0
            do i=0,dims
                do j=1,ns
                    do k=1,int(npx)
                        rl(j,k)=(x(i,j)-nx(k))/l
                    enddo
                enddo
                dl=rl**2
                i0=i*ns+1
                i1=i0+ns-1
                if (i.eq.0.and.gh.eq.0) then
        !            print *,'CV-Entro Matern'
                    call matern(dl,size(dl,1),size(dl,2))
                else
                    if (i.ge.1.and.gh.eq.0.or.i.eq.0.and.gh.eq.1) then
        !                print *,'covar vector calling deriv(1)'
                        call matderiv(1,size(dl,1),size(dl,2))
                        tmat=mat
                        mat=mat*(2*rl)/l
                    else
                        if(i.ge.1.and.gh.eq.1.or.i.eq.0.and.gh.eq.2) then
        !                    print *,'covar vector calling deriv(2) for Kriging Gradients'
                            call matderiv(2,size(dl,1),size(dl,2))
                            tmat2=mat
                            mat=mat*(-4*rl**2/l**2)+tmat*(-2/l**2)
                        else
                            if(i.ge.1.and.gh.eq.2) then
        !                        print *,'covar vector calling deriv(3)for Kriging Hessian'
                                call matderiv(1,size(dl,1),size(dl,2))
                                tmat=mat
                                call matderiv(3,size(dl,1),size(dl,2))
                                mat=mat*8*rl**4/l**3+tmat2*(16*rl**2/l**3 + 4/l**2)+tmat*(4/l**3)
                            endif
                        endif
                    endif
                endif
                cv(i0:i1,:)=mat
                !Write (6,*) cv
            enddo
        END
