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
            real*8 B(m_t,npx,nInter,nInter),A(m_t,m_t),tsum,ddottemp(npx),tcv(npx,m_t) !AF contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
            integer IPIV(m_t),INFO,i,j,iter,nInter,gh ! ipiv the pivot indices that define the permutation matrix
            A=full_R
            B=CV
            tsum=sum(rones(1:iter))
            CALL DGESV_(size(A,1), size(B,2),A,size(A,2),&
                    IPIV,B,size(B,1),INFO )
            do j=1,npx
                if (gh.eq.0) then
                    tcv=transpose(cv(:,:,1,1))
                    var=0
                    do i=1,m_t
                        var(j)=var(j)+B(i,j,1,1)*CV(i,j,1,1)
                    enddo
                    var=1-var
                    ddottemp(j)=dot_product(tcv(j,:),rones)
                    var=var+(1-ddottemp)**2/tsum
                    pred(j) = sb + dot_product(tcv(j,:),Kv)
                    sigma=1.96*sqrt(abs(var*variance))
!                   write(6,*) 'pred:',k,j,l,pred(j),var,variance, &
!                               sigma, lh
                else
                    if (gh.eq.1) then
                        do k=1,nInter
                            tcv=transpose(cv(:,:,k,1))
                            gpred(j,k) = dot_product(tcv(j,:),Kv)
                            sigma=1.96*sqrt(2*abs(var*variance))
                            ! write(6,*) 'pred Grad:',k,j,l,gpred(j,k), &
                            !     var,variance,sigma, lh,tcv
                        enddo
!                       write(6,*) 'final cv', cv(:,:,:,1)
!                       write(6,*) 'final Kv',kv
!                       write (6,*) 'pred grad:',gpred
                    else
                        do k=1,nInter
                            do i=1,nInter
                                tcv=transpose(cv(:,:,k,i))
                                hpred(j,k,i) = dot_product(tcv(j,:),Kv)
                                sigma=1.96*sqrt(2*abs(var*variance))
                                ! write(6,*) 'pred Hess:',k,j,l,hpred(j,k), &
                                !     var,variance,sigma, lh, tcv
                            enddo
                        enddo
                    endif
                endif
            enddo
        END
