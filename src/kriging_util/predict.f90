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
            real*8 B(m_t,npx,nInter),A(m_t,m_t),tsum,ddottemp(npx),tcv(npx,m_t) !AF contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
            integer IPIV(m_t),INFO,i,j,iter,nInter,gh ! ipiv the pivot indices that define the permutation matrix
            ! write(6,*) 'Predict: '
            ! variance=dot_product(Ys,Kv)/m_t
            ! if (detr<0) then
            !     lh=-variance*exp(log(-detr)/m_t)
            ! else
            !     lh=variance*exp(log(detr)/m_t)
            ! endif
            A=full_R
            B=CV
            tsum=sum(rones(1:iter))
            CALL DGESV_(size(A,1), size(B,2),A,size(A,2),&
                    IPIV,B,size(B,1),INFO )
            do k=1,nInter
                tcv=transpose(cv(:,:,k))
                do j=1,npx
                    var=0
                    do i=1,m_t
                        var(j)=var(j)+B(i,j,k)*CV(i,j,k)
                    enddo
                    var=1-var
                    ddottemp(j)=dot_product(tcv(j,:),rones)
                    var=var+(1-ddottemp)**2/tsum
                    if (gh.eq.0) then
                        pred(j) = sb + dot_product(tcv(j,:),Kv)
                        sigma=1.96*sqrt(abs(var*variance))
!                       write(6,*) 'pred:',k,j,l,pred(j),var,variance, &
!                                   sigma, lh
                    else
                        pred(j) = dot_product(tcv(j,:),Kv)
                        sigma=1.96*sqrt(2*abs(var*variance))
!                       write(6,*) 'pred Grad:',k,j,l,pred(j),var,variance, &
!                                   sigma, lh
                    endif
                enddo
            enddo

            ! do k=1,nInter
            !     do j=1,npx
            !         do i=1,m_t
            !             var(j)=var(j)+B(i,j,k)*CV(i,j,k)
            !         enddo
            !     enddo
            ! enddo
            ! do i=1,npx
            !     ddottemp(i)=dot_product(tcv(i,:),rones)
            ! enddo
            !write(6,*) "Kv: ",Kv
            !write(6,*) "ys: ",Ys
!           write(6,*) 'l:',l
!           write(6,*) 'var: ',var
!           write(6,*) 'variance: ',variance
!           write(6,*) 'sigma: ',sigma
            !write(6,*) 'll: ',ll
        END
