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
        SUBROUTINE covarMatrix()
            use globvar
            integer i,j,i0,i1,j0,j1,k,kl
            real*8 tmat(ns,ns),diffx(ns,ns),diffx0(ns,ns),c(ns,ns)
            deallocate (rl,dl,mat)
            allocate (rl(NS,NS),dl(NS,NS),mat(NS,NS))
            Write (6,*) 'Covar Matrix'
            full_R=0
            c=exp(-sqrt((2.0*p+1.0)*dl))
            tmat=0
            ! Covariant Matrix in kriging
            do i=1,dims
                do k=1,NS
                    do kl=1,NS
                        rl(k,kl)=(x(i,k)-x(i,kl))/l
                    end do
                end do
                dl=rl**2
                call matern(dl,ns,ns)
                tmat=tmat+mat
            end do
            full_R(1:NS,1:NS)=tmat+iden*eps
            ! Covariant matrix in Gradient Enhanced Kriging
            do i=1,dims
                do k=1,NS
                    do kl=1,NS
                        rl(k,kl)=(x(i,k)-x(i,kl))/l
                    end do
                end do
                dl=rl**2
                diffx=2.0*rl/l
                diffx0=-2.0*rl/l
                i0=i*ns+1
                i1=i0+ns-1
                do j=i,dims+1
                    j0=(j-1)*ns+1
                    j1=j0+ns-1
                    if(j.eq.1) then
                        call matderiv(1,ns,ns)
                        tmat=mat
                        mat=mat*diffx
                    else
                        call matderiv(2,ns,ns)
                        mat=mat*diffx*diffx0-tmat*(2/l**2)
                    endif
                    full_R(i0:i1,j0:j1)=mat
                    if (j.ne.dims+1) then
                        full_R(j0:j1,i0:i1)=transpose(mat)
                    else
                        full_R(i0:i1,j0:j1)=full_R(i0:i1,j0:j1)+iden*eps
                    endif
                enddo
            enddo
            write (6,*) 'full_R: ',full_R
            Write (6,*) 'full_R size: ',size(full_R)
            Write (6,*) 'full_R shape: ',shape(full_R)
        END