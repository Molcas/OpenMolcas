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
    integer i,j,i0,i1,j0,j1
    real*8 tmat(ns,ns),diffx(ns,ns),diffx0(ns,ns),c(ns,ns)
    deallocate (rl,dl,mat)
    allocate (rl(NS,NS),dl(NS,NS),mat(NS,NS))
    do i=1,NS
        do j=1,NS
            rl(i,j)=(x(i)-x(j))/l
        end do
    end do
    dl=rl**2
    diffx=2.0*rl/l
    diffx0=-2.0*rl/l
    der=0
    full_R=0
    c=exp(-sqrt((2.0*p+1.0)*dl))
    do i=0,dims
        i0=i*ns+1
        i1=i0+ns-1
        if (i.eq.0) then
            der(1,:)=0
        else
            der(1,:)=[0,1]
        endif
        do j=i,dims
            j0=j*ns+1
            j1=j0+ns-1
            if(i.eq.0.and.j.eq.0) then
                call matern(dl,ns,ns)
            else
                if(i.eq.0.and.j.eq.1) then
                    call matderiv(1,ns,ns)
                    tmat=mat
                    mat=mat*diffx
                else
                    if(i.eq.1.and.j.eq.1) then
                        call matderiv(2,ns,ns)
                        mat=mat*diffx*diffx0-tmat*(2/l**2)
                    endif
                endif
            endif
!            print *,'CovarMatrix mat'
!            call printsmat(mat,size(mat,1),size(mat,2))
            full_R(i0:i1,j0:j1)=transpose(mat)
            if (i.ne.j) then
                full_R(j0:j1,i0:i1)=mat
            else
                full_R(i0:i1,j0:j1)=full_R(i0:i1,j0:j1)+iden*eps
            endif
        enddo
    enddo
END
