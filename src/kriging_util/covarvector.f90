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
        SUBROUTINE covarVector(gh,iter,nInter)
            use globvar
#include "stdalloc.fh"
            integer i,i0,i1,j,gh,iter,nInter
            real*8 sdiffx,sdiffx0,sdiffxk
            real*8, Allocatable ::  diffx(:,:),diffx0(:,:), diffxk(:,:)
!
            Call mma_Allocate(diffx,iter,npx,label="diffx")
            Call mma_Allocate(diffx0,iter,npx,label="diffx0")
            Call mma_Allocate(diffxk,iter,npx,label="diffxk")
!
            cv = 0
            i0 = 0
            call defdlrl(iter,nInter)
!           write(6,*) 'x: ',x
!           write(6,*) 'nx: ',nx
! Covariant Vector in kriging - First part of eq (4) in ref.
!
            if (gh.eq.0) then
!
                call matern(dl, cv(1:iter,:,1,1), iter, npx)
                call matderiv(1, dl, cvMatFDer, iter, npx)
                do i=1,nInter
!       1st derivatives second part of eq. (4)
                    diffx(:,:) = 2.0D0*rl(:,:,i)/l(i)
                    i0 = i*iter + 1
                    i1 = i0 + iter - 1
                    cv(i0:i1,:,1,1) = cvMatFder * diffx
                enddo
                ! write (6,*) 'CV-mat',cv
! Covariant vector in Gradient Enhanced Kriging
!
            else if(gh.ge.1) then
!
                ! print *,'covar vector calling deriv(2) for Kriging Gradients'
                call matderiv(1, dl, cvMatFder, iter, npx)
                ! Call RecPrt('cvMatFder',' ',cvMatFder,iter,npx)
                call matderiv(2, dl, cvMatSder, iter, npx)
                ! Call RecPrt('cvMatSder',' ',cvMatSder,iter,npx)
                do i=1,nInter
                    diffx(:,:) = 2.0D0*rl(:,:,i)/l(i)
                    cv(1:iter,:,i,1) = -cvMatFder * diffx
                    do j = 1,nInter
                        ! if (j.eq.1) cv(1:iter,:,i,1) = -cvMatFder * diffx
                        j0 = j*iter + 1
                        j1 = j0+iter - 1
                        diffx0(:,:) = -2.0D0*rl(:,:,j)/l(j)
                        if (i.eq.j) Then
                           cv(j0:j1,:,i,1) = cvMatSder * diffx*diffx0 - cvMatFder*(2/(l(i)*l(j)))
                        else
                           cv(j0:j1,:,i,1) = cvMatSder * diffx*diffx0
                        end if
                        ! Write (6,*) 'cvMatSder - Krig Grad: ',cvMatSder
                        ! write (6,*) 'CV',cv(:,:,:,1)
                    enddo
                enddo
                ! Write (6,*) 'CV - Krig Grad: ',cv
!
            else if(gh.eq.2) then
!
!                    print *,'covar vector calling deriv(3) for Kriging Hessian'
                ! anAI = .False.
                call matderiv(1, dl, cvMatFder, iter, npx)
                call matderiv(2, dl, cvMatSder, iter, npx)
                call matderiv(3, dl, cvMatTder, iter, npx)
                do i = 1, nInter
                    diffx(:,:) = 2.0D0*rl(:,:,i)/l(i)
                    sdiffx = 2.0D0/l(i)**2
                    do j = 1, nInter
                        diffx0(:,:) = -2.0D0*rl(:,:,j)/l(j)
                        sdiffx0 = 2.0D0/l(j)**2
                        if (i.eq.j) Then
                           cv(1:iter,:,i,j) = cvMatSder * diffx*diffx0 - cvMatFder*2.0D0/(l(i)*l(j))
                        else
                           cv(1:iter,:,i,j) = cvMatSder * diffx*diffx0
                        end if
                        do k = 1, nInter
                            diffxk(:,:) = - 2.0D0*rl(:,:,k)/l(k)
                            sdiffxk = 2.0D0/l(i)**2
                            k0 = k*iter + 1
                            k1 = k0+iter - 1
                            if (i.eq.j.and.j.eq.k) then
                                cv(k0:k1,:,i,j) = (cvMatTder*diffx0**3 + 3.0D0*cvMatSder*diffx0*sdiffx0)
                            else if (i.eq.j) then
                                cv(k0:k1,:,i,j) = cvMatTder*diffx0**3 + cvMatSder*diffx0*sdiffx
                            else if (i.eq.k) then
                                cv(k0:k1,:,i,j) = cvMatTder*diffxk**3 + cvMatSder*diffxk*sdiffx
                            else if (j.eq.k) then
                                cv(k0:k1,:,i,j) = cvMatTder*diffx0**3 + cvMatSder*diffxk*sdiffxk
                            else
                                cv(k0:k1,:,i,j) = cvMatTder*diffx*diffx0*diffxk
                            endif
                        enddo
                    enddo
                enddo
                !Write (6,*) 'CV - Krig Hessian: ',cv
            else
                Write (6,*) ' Illegal value of gh:',gh
                Call Abend()
            endif
            ! Write (6,*) 'CV shape: ',shape(CV)
            ! write (6,*) 'CV: ',CV
!
            Call mma_deallocate(diffx)
            Call mma_deallocate(diffx0)
            Call mma_deallocate(diffxk)
!
        END Subroutine covarvector
!
        SUBROUTINE defdlrl(iter,nInter)
            use globvar
            integer i,j,iter,nInter
            dl=0
            do i=1,nInter
                do j=1,iter
                    do k=1,int(npx)
                        rl(j,k,i) = (x(i,j) - nx(i,k))/l(i)
                        ! write (6,*) i,j,k,'rl,x,nx',rl(j,k,i),x(i,j),nx(i,k),l(i)
                    enddo
                enddo
                !write(6,*) 'CV-rl',i,rl
                dl(:,:) = dl(:,:) + rl(:,:,i)**2
            enddo
            ! write (6,*) 'rl',rl
            !isdefdlrl = .True.
        END Subroutine defdlrl
