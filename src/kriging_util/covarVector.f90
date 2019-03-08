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
            do i=1,dims
                do j=1,ns
                    do k=1,int(npx)
                        rl(j,k)=(x(i,j)-nx(k))/l
                    enddo
                enddo
                dl=rl**2
                i0=(i-1)*ns+1
                i1=i0+ns-1
                Write (6,*) 'Covar Vector 39', i, io, i1
                if (i.eq.1.and.gh.eq.0) then
        !            print *,'CV-Entro Matern'
                    call matern(dl,size(dl,1),size(dl,2))
                    Write (6,*) 'Covar Vector 43'
                else
                    if (i.ge.2.and.gh.eq.0.or.i.eq.1.and.gh.eq.1) then
        !                print *,'covar vector calling deriv(1)'
                        Write (6,*) 'Covar Vector 50'
                        call matderiv(1,size(dl,1),size(dl,2))
                        Write (6,*) 'Covar Vector 50'
                        tmat=mat
                        Write (6,*) 'Covar Vector 50'
                        mat=mat*(2*rl)/l
                        Write (6,*) 'Covar Vector 50'
                    else
                        if(i.ge.2.and.gh.eq.1.or.i.eq.1.and.gh.eq.2) then
        !                    print *,'covar vector calling deriv(2) for Kriging Gradients'
                            call matderiv(2,size(dl,1),size(dl,2))
                            tmat2=mat
                            mat=mat*(-4*rl**2/l**2)+tmat*(-2/l**2)
                            Write (6,*) 'Covar Vector 57'
                        else
                            if(i.ge.2.and.gh.eq.2) then
        !                        print *,'covar vector calling deriv(3)for Kriging Hessian'
                                call matderiv(1,size(dl,1),size(dl,2))
                                tmat=mat
                                call matderiv(3,size(dl,1),size(dl,2))
                                mat=mat*8*rl**4/l**3+tmat2*(16*rl**2/l**3 + 4/l**2)+tmat*(4/l**3)
                                Write (6,*) 'Covar Vector 65'
                            endif
                        endif
                    endif
                endif
                cv(i0:i1,:)=mat
                !Write (6,*) cv
            enddo
        END
