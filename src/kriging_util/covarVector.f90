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
            integer i,i0,i1,j,gh,k,iter,nInter
            ! real*8 tmat(iter,npx),tmat2(iter,npx), &
            real*8 m(iter,npx),diffx(iter,npx),diffx0(iter,npx), &
                   dl(iter,npx)
            ! deallocate (dl,mat)
            ! allocate (dl(iter,npx),mat(iter,npx))
            cv = 0
!           write(6,*) 'x: ',x
!           write(6,*) 'nx: ',nx
! Covariant Vector in kriging - First part of eq (4) in ref.
            if (gh.eq.0) then
                ! tmat=0 ! to be removed
                ! tmat2=0
                dl=0
                do i=1,nInter
                    do j=1,iter
                        do k=1,int(npx)
                            rl(j,k,i) = (x(i,j) - nx(i,k))/l(i)
                        enddo
                    enddo
                    !write(6,*) 'CV-rl',i,rl
                    dl = dl + rl(:,:,i)**2
                    !write(6,*) 'CV-Mat',i,mat
                    !tmat=tmat+mat
                enddo
!               write(6,*) 'CV-rl',rl
!               write(6,*) 'CV-d',dl
                call matern(dl, m, iter, npx)
                cv(1:iter,:,1,1) = m
!               write (6,*) 'cv-gh,m',gh,m
                call matderiv(1, dl, m, iter, npx)
                cvMatFder = m
                call matderiv(2, dl, m, iter, npx)
                cvMatSder = m
                call matderiv(3, dl, m, iter, npx)
                cvMatTder = m
                ! write (6,*) 'CV-mat',cv
            endif
! Covariant vector in Gradient Enhanced Kriging
            !do i=1,nInter
                ! do j=1,iter
                !     do k=1,int(npx)
                !         rl(j,k)=(x(i,j)-nx(i,k))/l(i)
                !     enddo
                ! enddo
                ! dl=rl**2
                ! if (gh.eq.0) then
                !     call matern(dl,size(dl,1),size(dl,2))
                !     cv(1:iter,:,i)=mat
                ! endif
                ! write(6,*) 'diffx',diffx
                ! write(6,*) 'diffx0',diffx0
                if (gh.eq.0) then
                    do i=1,nInter
!       1st derivatives second part of eq. (4)
                        diffx = 2.0*rl(:,:,i)/l(i)
                        i0 = i*iter + 1
                        i1 = i0 + iter - 1
!                       write (6,*) 'i,i0,i1',i,i0,i1
                        !call matderiv(1, dl, m, iter, npx)
                        !tmat = m
                        m = cvMatFder * diffx
                        cv(i0:i1,:,1,1) = m
                        !write(6,*) 'tmat',tmat
!                       write(6,*) 'gh,m',gh,m
!                       write(6,*) 'cv',cv(:,:,1,1)
                    enddo
                ! else
                !     i0 = (i-1)*iter + 1
                !     i1 = i0+iter - 1
                endif
                !write(6,*) 'tmat,m',tmat,m
                !Write (6,*) 'CV - Deriv: ',mat
                if(gh.ge.1) then
!                    print *,'covar vector calling deriv(2) for Kriging Gradients'
                    ! if (gh.eq.1) then
                    !     write(6,*) 'm',m
                    !     cv(i0:i1,:,i,1) = -m
                    !     i0 = i*iter + 1
                    !     i1 = i0 + iter - 1
                    ! endif
                    ! write(6,*) 'dl',dl
                    ! call matderiv(2, dl, m, iter, npx)
                    ! tmat2 = m
                    ! write(6,*) 'tmat2',tmat2
                    do i=1,nInter
                        diffx = 2.0*rl(:,:,i)/l(i)
                        do j = 1,nInter
                            if (j.eq.1) cv(1:iter,:,i,1) = -cvMatFder * diffx
                            j0 = j*iter + 1
                            j1 = j0+iter - 1
!                           write(6,*) 'i,j',i,j
                            diffx0 = -2.0*rl(:,:,j)/l(j)
!                           write(6,*) 'diffx',diffx
!                           write(6,*) 'diffx0',diffx0
                            m = cvMatSder * diffx*diffx0
                            if (i.eq.j) m = m - cvMatFder*(2/(l(i)*l(j)))
                            cv(j0:j1,:,i,1) = m
!                           Write (6,*) 'm - Krig Grad: ',m
!                           write (6,*) 'CV',cv(:,:,:,1)
                        enddo
                    enddo
                    !Write (6,*) 'CV - Krig Grad: ',cv
                endif
                if(gh.eq.2) then
!                    print *,'covar vector calling deriv(3)for Kriging Hessian'
                    cv(i0:i1,:,i,i) = m
                    i0 = i*iter + 1
                    i1 = i0+iter - 1
                    call matderiv(3, dl, m, iter, npx)
                    m = cvMatTder * diffx*diffx0**2 + cvMatSder*(4*diffx**2/l(i) + &
                        4/l(i)**2) + cvMatFder*(4/l(i)**3)
                    !Write (6,*) 'CV - Krig Hessian: ',mat
                endif
            ! enddo
            ! Write (6,*) 'CV shape: ',shape(CV)
            ! write (6,*) 'CV: ',CV
        END
