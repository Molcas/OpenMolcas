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
            integer i,i0,i1,j,gh,iter,nInter
            ! real*8 tmat(iter,npx),tmat2(iter,npx), &
            real*8 m(iter,npx),diffx(iter,npx),diffx0(iter,npx), &
                diffxk(iter,npx),sdiffx,sdiffx0!, & dl(iter,npx)
            ! deallocate (dl,mat)
            ! allocate (dl(iter,npx),mat(iter,npx))
            cv = 0
            i0 = 0
!           write(6,*) 'x: ',x
!           write(6,*) 'nx: ',nx
! Covariant Vector in kriging - First part of eq (4) in ref.
            if (gh.eq.0) then
                ! tmat=0 ! to be removed
                ! tmat2=0
!               write(6,*) 'CV-rl',rl
!               write(6,*) 'CV-d',dl
                call defdlrl(iter,nInter)
                call matern(dl, m, iter, npx)
                cv(1:iter,:,1,1) = m
!               write (6,*) 'cv-gh,m',gh,m
                call matderiv(1, dl, m, iter, npx)
                cvMatFder = m
                ! call matderiv(2, dl, m, iter, npx)
                ! cvMatSder = m
                ! call matderiv(3, dl, m, iter, npx)
                ! cvMatTder = m
                do i=1,nInter
!       1st derivatives second part of eq. (4)
                    diffx = 2.0*rl(:,:,i)/l(i)
                    i0 = i*iter + 1
                    i1 = i0 + iter - 1
                    m = cvMatFder * diffx
                    cv(i0:i1,:,1,1) = m
                enddo
                ! write (6,*) 'CV-mat',cv
            endif
! Covariant vector in Gradient Enhanced Kriging
            if(gh.ge.1) then
!                    print *,'covar vector calling deriv(2) for Kriging Gradients'
                !if (.NOT. isdefdlrl) then
                call defdlrl(iter,nInter)
                call matderiv(1, dl, m, iter, npx)
                cvMatFder = m
                call matderiv(2, dl, m, iter, npx)
                cvMatSder = m
                !endif
                do i=1,nInter
                    diffx = 2.0*rl(:,:,i)/l(i)
                    cv(1:iter,:,i,1) = -cvMatFder * diffx
                    do j = 1,nInter
                        ! if (j.eq.1) cv(1:iter,:,i,1) = -cvMatFder * diffx
                        j0 = j*iter + 1
                        j1 = j0+iter - 1
!                           write(6,*) 'i,j',i,j
                        diffx0 = -2.0*rl(:,:,j)/l(j)
!                           write(6,*) 'diffx',diffx
!                           write(6,*) 'diffx0',diffx0
                        m = cvMatSder * diffx*diffx0
                        if (i.eq.j) m = m - cvMatFder*(2/(l(i)*l(j)))
                        cv(j0:j1,:,i,1) = m
                        ! Write (6,*) 'cvMatSder - Krig Grad: ',cvMatSder
!                           write (6,*) 'CV',cv(:,:,:,1)
                    enddo
                enddo
                !Write (6,*) 'CV - Krig Grad: ',cv
            endif
            if(gh.eq.2) then
!                    print *,'covar vector calling deriv(3) for Kriging Hessian'
                call defdlrl(iter,nInter)
                call matderiv(1, dl, m, iter, npx)
                cvMatFder = m
                call matderiv(2, dl, m, iter, npx)
                cvMatSder = m
                Write (6,*) 'cvMatSder - Krig Grad(hess): ',cvMatSder
                call matderiv(3, dl, m, iter, npx)
                cvMatTder = m
                ! write (6,*) 'dl',dl
                ! write (6,*) '3th der',cvMatTder
                do i = 1, nInter
                    diffx = 2.0*rl(:,:,i)/l(i)
                    sdiffx = 2.0/l(i)**2
                    do j = 1, nInter
                        diffx0 = -2.0*rl(:,:,j)/l(j)
                        sdiffx0 = 2/l(j)**2
                        m = cvMatSder * diffx*diffx0
                        if (i.eq.j) m = m - cvMatFder*2/(l(i)*l(j))
                        cv(1:iter,:,i,j) = m
                        do k = 1, nInter
                            diffxk = -2.0*rl(:,:,k)/l(k)
                            k0 = k*iter + 1
                            k1 = k0+iter - 1
                            if (i.eq.j.and.j.eq.k) then
                                m = cvMatTder*diffx0**3 + 3*cvMatSder*diffx0*sdiffx0
                                ! write(6,*) 'i=j=k',i,j,k
                                ! write(6,*) 'cvMatTder*diffx0**3',cvMatTder*diffx0**3
                            else
                                if (i.eq.j) then
                                    m = cvMatTder*diffx0**3 + cvMatSder*diffx0*sdiffx
                                    !write(6,*) 'i=j!=k',i,j,k
                                else
                                    if (i.eq.k) then
                                        m = cvMatTder*diffxk**3 + cvMatSder*diffxk*sdiffx
                                        !write(6,*) 'i=K!=J',i,j,k
                                    else
                                        if (j.eq.k) then
                                            m = cvMatTder*diffx0**3 + cvMatSder*diffxk*sdiffx0
                                            !write(6,*) 'i=j!=k',i,j,k
                                        else
                                            m = -cvMatTder*diffx*diffx0*diffxk
                                            ! write(6,*) 'i!=j!=k',i,j,k
                                            ! write (6,*) m
                                        endif
                                    endif
                                endif
                            endif
                            !write(6,*) 'm',m
                            ! m = cvMatTder * diffx*diffx0**2 + cvMatSder*(4*diffx**2/l(i) + &
                            !     4/l(i)**2) + cvMatFder*(4/l(i)**3)
                            cv(k0:k1,:,i,j) = m
                        enddo
                    enddo
                enddo
                !Write (6,*) 'CV - Krig Hessian: ',cv
            endif
            ! Write (6,*) 'CV shape: ',shape(CV)
            ! write (6,*) 'CV: ',CV
        END
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
                    dl = dl + rl(:,:,i)**2
                enddo
                ! write (6,*) 'rl',rl
                !isdefdlrl = .True.
        END