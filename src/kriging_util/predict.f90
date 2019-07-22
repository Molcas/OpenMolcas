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
        SUBROUTINE predict(gh,iter,nInter)
            use globvar
            real*8 B(m_t),tsum!,ddottemp(npx)!,tcv(npx,m_t)
            real*8 A(m_t,m_t)!AF contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
            integer IPIV(m_t),INFO
            integer i,j,iter,nInter,gh ! ipiv the pivot indices that define the permutation matrix
!
            do j=1,npx
                if (gh.eq.0) then
! ---------------calculations for  dispersion ----------------------------
                    ! ----------------Old calculations --P1
                    A = full_R
                    B = CV(:,j,1,1)
                    CALL DGESV_(m_t, 1,A,m_t,IPIV,B,m_t,INFO )
                    !-----------------New
                    ! Call RecPrt('full_Rinv',  ' ',full_Rinv,m_t,m_t)
                    ! B = CV(:,j,1,1)
                    ! B = matmul(B,full_Rinv)
                    !--------------------
                    ! write (6,*) 'R^(-1)*CV', B
                    ! write (6,*) 'size', size(B)
                    var(j) = 1 - dot_product(B,CV(:,j,1,1))
                    tsum = sum(rones(1:iter))
                    B = cv(:,j,1,1)
                    var(j)=var(j)+(1-dot_product(B,rones))**2/tsum
                    sigma(j)=1.96*sqrt(abs(var(j)*variance))
! -------------------------------
! -----Predictiong the Energy
                    pred(j) = sb + dot_product(B,Kv)
                !   write(6,*) 'pred:',pred(j)
                !   write(6,*) 'var:',var
                !   write(6,*) 'variance',variance
                !   write(6,*) 'sigma',sigma,'lh',lh
                !   write(6,*) 'tcv(j,:)',B
                !   write(6,*) 'Kv',Kv
                else
                    if (gh.eq.1) then
                        ! sigma(j)=1.96*sqrt(2*abs(var*variance))
                        do k=1,nInter
                            B = cv(:,j,k,1)
                            gpred(j,k) = dot_product(B,Kv)
                            ! write(6,*) 'pred Grad:',k,j,l,gpred(j,k), &
                            !     var,variance,sigma, lh,tcv
                        enddo
                    !   write(6,*) 'final cv', cv(:,:,:,1)
                    !   write(6,*) 'final Kv',kv
                    !   write (6,*) 'pred grad(gpred):',gpred
                    else
                        ! sigma(j)=1.96*sqrt(2*abs(var*variance))
                        ! write(6,*) 'kv: ',kv
! Predicting the Hessian gh = 2
                        do k=1,nInter
                            do i=1,nInter
                                B = cv(:,j,i,k)
                                ! write(6,*) 'tcv', i,k,tcv
                                !Call RecPrt('Update_: tcv',' ',tcv,npx,m_t)
                                hpred(j,k,i) = dot_product(B, Kv)
                                !write (6,*) 'partial hpred',hpred(j,k,i)
                                ! write(6,*) 'pred Hess:',k,j,l,hpred(j,k), &
                                !     var,variance,sigma, lh, tcv
                            enddo
                        enddo
                        ! write (6,*) 'pred hess(hpred):',hpred
                    endif
                endif
            enddo
        END
