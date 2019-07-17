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
            real*8 diffx(iter,iter),diffx0(iter,iter), &
                    matFder(iter,iter),matSder(iter,iter),r(iter,iter,nInter), &
                    d(iter,iter),m(iter,iter),iden(iter,iter)!,c(iter,iter)
            full_R = 0
            d = 0
            diffx = 0
            diffx0 = 0
            call miden(iden,iter)
! Covariant Matrix in kriging
            do i=1,nInter
                do k=1,iter
                    do kl=1,iter
                        r(k,kl,i)=(x(i,k)-x(i,kl))/l(i)
                    end do
                end do
                d = d + r(:,:,i)**2
            end do
    !Matern Function
            Call matern(d, m, iter, iter)
! Writing the covariant matrix in GEK (eq 2 of DOI 10.1007/s00366-015-0397)
            full_R(1:iter,1:iter) = m + iden*eps
    !Matern first derivative
            call matderiv(1, d, m, iter, iter)
            matFder = m
    !Matern second derivative
            call matderiv(2, d, m, iter, iter)
            matSder = m
! Covariant matrix in Gradient Enhanced Kriging (eq 2 of DOI 10.1007/s00366-015-0397)):
!
    ! First line and first column derivative in Psi matrix
            do i=1,nInter
                i0=i*iter+1
                i1=i0+iter-1
                diffx0 = -2.0*r(:,:,i)/l(i)
                m = matFder*diffx0
    !  Writing the 1st row of 1st derivatives with respect the coordinates
                full_R(1:iter,i0:i1) = m
    !  Writing the column of derivatives
                full_R(i0:i1,1:iter) = transpose(m)
            enddo
    ! Second derivatives
            do i = 1,nInter
                i0 = i*iter+1
                i1 = i0+iter-1
                do j = i,nInter
                    j0 = j*iter+1
                    j1 = j0+iter-1
                    diffx = 2.0*r(:,:,j)/l(j)
                    diffx0 = -2.0*r(:,:,i)/l(i)
                    m = matSder*diffx*diffx0
    !   if differentiating twice on the same dimension
                    if (i.eq.j) m = m - matfder*(2/(l(i)*l(j)))
    !   Writing the second derivatives in eq(2)
                    full_R(i0:i1,j0:j1) = m
                    if (i.ne.j) then
                        full_R(j0:j1,i0:i1) = transpose(m)
                    else
                        full_R(i0:i1,j0:j1) = full_R(i0:i1,j0:j1) + iden*eps2
                    endif
                enddo
            enddo
!           definig full_r has srictly possitive define sec. 3 of
!           DOI: 10.1615/Int.J.UncertaintyQuantification.2013006809
            ! full_R = abs(full_R)
        !   Call RecPrt('full_r',  ' ',full_R,m_t,m_t)
        END