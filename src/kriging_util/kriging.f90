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
        SUBROUTINE predict()
            use globvar
            real*8 B(m_t,npx),A(m_t,m_t),tsum,ddottemp(npx),tcv(npx,m_t) !AF contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
            integer IPIV(m_t),INFO,i,j ! ipiv the pivot indices that define the permutation matrix
            write(6,*) 'Predict'
            variance=dot_product(Ys,Kv)/m_t
            A=full_R
            B=CV
            CALL DGESV_(size(A,1), size(B,2),A,size(A,2),&
                    IPIV,B,size(B,1),INFO )
            var=0
            do j=1,npx
                do i=1,m_t
                    var(j)=var(j)+B(i,j)*CV(i,j)
                enddo
            enddo
            var=1-var
            tcv=transpose(cv)
            tsum=sum(rones(1:ns))
            do i=1,npx
                ddottemp(i)=dot_product(tcv(i,:),rones)
            enddo
            var=var+(1-ddottemp)**2/tsum
            do i=1,npx
                pred(i) = sb + dot_product(tcv(i,:),Kv)
                write(6,*) 'pred', i, pred(i)
            enddo
            sigma=1.96*sqrt(var*variance)
            if (detr<0) then
                ll=-variance*exp(log(-detr)/m_t)
            else
                ll=variance*exp(log(detr)/m_t)
            endif
            write(6,*) 'l ',l
            write(6,*) 'variance ',variance
            write(6,*) 'll ',ll
        END