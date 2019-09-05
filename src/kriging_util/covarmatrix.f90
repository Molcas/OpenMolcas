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
#include "stdalloc.fh"
            integer i,j,i0,i1,j0,j1,k,kl,iter,nInter
            Real*8, Allocatable :: diffx(:,:), diffx0(:,:), matFder(:,:),&
                                   matSder(:,:), r(:,:,:), d(:,:)
!
            Call mma_Allocate(diffx,iter,iter,Label="diffx")
            Call mma_Allocate(diffx0,iter,iter,Label="diffx0")
            Call mma_Allocate(matFder,iter,iter,Label="matFder")
            Call mma_Allocate(matSder,iter,iter,Label="matSder")
            Call mma_Allocate(r,iter,iter,nInter,Label="r")
            Call mma_Allocate(d,iter,iter,Label="d")
!
            full_R = 0
            d = 0
            diffx = 0
            diffx0 = 0
!
! Covariant Matrix in kriging
            do i=1,nInter
                do k=1,iter
                    do kl=1,iter
                        r(k,kl,i)=(x(i,k)-x(i,kl))/l(i)
                    end do
                end do
                d(:,:) = d(:,:) + r(:,:,i)**2
            end do
!
    !Matern Function
            Call matern     (d, full_R(1:iter,1:iter), iter, iter)
!
! Writing the covariant matrix in GEK (eq 2 of DOI 10.1007/s00366-015-0397)
!
    !Matern first derivative
            call matderiv(1, d, MatFder, iter, iter)
! Covariant matrix in Gradient Enhanced Kriging (eq 2 of DOI 10.1007/s00366-015-0397)):
!
    ! First line and first column derivative in Psi matrix
            do i=1,nInter
                i0=i*iter+1
                i1=i0+iter-1
                diffx0(:,:) = -2.0D0*r(:,:,i)/l(i)
    !  Writing the 1st row of 1st derivatives with respect the coordinates
                full_R(1:iter,i0:i1) = matFDer*diffx0
    !  Writing the column of derivatives
                full_R(i0:i1,1:iter) = transpose(full_R(1:iter,i0:i1))
            enddo
!
    !Matern second derivative
            call matderiv(2, d, matSder, iter, iter)
!
    ! Second derivatives
            do i = 1,nInter
                i0 = i*iter+1
                i1 = i0+iter-1
                do j = i,nInter
                    j0 = j*iter+1
                    j1 = j0+iter-1
                    diffx(:,:)  =  2.0D0*r(:,:,j)/l(j)
                    diffx0(:,:) = -2.0D0*r(:,:,i)/l(i)
    !   if differentiating twice on the same dimension
                    if (i.eq.j) Then
                       full_R(i0:i1,j0:j1) = matSder*diffx*diffx0 - matfder*(2/(l(i)*l(j)))
                    else
                       full_R(i0:i1,j0:j1) = matSder*diffx*diffx0
                    end if
    !   Writing the second derivatives in eq(2)
                    if (i.ne.j) full_R(j0:j1,i0:i1) = transpose(Full_r(i0:i1,j0:j1))
                enddo
            enddo
!
!           Add constants to reflect the error in the energy and the
!           gradient, respectively.
!
            forall (j=1:iter) Full_R(j,j) = Full_R(j,j) + eps
            forall (j=iter+1:m_t) Full_R(j,j) = Full_R(j,j) + eps2
!
!           definig full_r has srictly possitive define sec. 3 of
!           DOI: 10.1615/Int.J.UncertaintyQuantification.2013006809
            ! full_R = abs(full_R)
        !   Call RecPrt('full_r',  ' ',full_R,m_t,m_t)
!
            Call mma_deallocate(diffx)
            Call mma_deallocate(diffx0)
            Call mma_deallocate(matFder)
            Call mma_deallocate(matSder)
            Call mma_deallocate(r)
            Call mma_deallocate(d)
        END
