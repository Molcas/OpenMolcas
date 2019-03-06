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
            full_R(1:NS,1:NS)=tmat
            do i=0,dims
                do k=1,NS
                    do kl=1,NS
                        rl(k,kl)=(x(i+1,k)-x(i+1,kl))/l
                    end do
                end do
                dl=rl**2
                diffx=2.0*rl/l
                diffx0=-2.0*rl/l
                i0=i*ns+1
                i1=i0+ns-1
                do j=i,dims
                    if (j.ge.1) then
                        j0=j*ns+1
                        j1=j0+ns-1
                        if(i.eq.0) then
                            call matderiv(1,ns,ns)
                            tmat=mat
                            mat=mat*diffx
                        else
                            call matderiv(2,ns,ns)
                            mat=mat*diffx*diffx0-tmat*(2/l**2)
                        endif
                        full_R(i0:i1,j0:j1)=transpose(mat)
                        if (i.ne.j) then
                            full_R(j0:j1,i0:i1)=mat
                        else
                            full_R(i0:i1,j0:j1)=full_R(i0:i1,j0:j1)+iden*eps
                        endif
                    endif
                enddo
            enddo
            !Write (6,*) full_R
        END