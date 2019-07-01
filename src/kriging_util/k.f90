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
        SUBROUTINE k(iter)
            use globvar
            real*8 B(m_t),A(m_t,m_t),diagA(m_t), iden(m_t,m_t)
            integer IPIV(m_t),INFO,iter,sign ! ipiv the pivot indices that define the permutation matrix
!
! Initiate B according to Eq. (6) of ref.
            B=0.0D0
            B(1:iter)=1.0D0
!
            call miden(iden,m_t)
            full_Rinv = iden
! Initiate A according to Eq. (2) of ref.
            A = full_r !in, coefficent matrix A, out factors L and U from factorization A=PLU on AX=B
!
            ! ----------------Old calculations --K1
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            rones=B
            !-----------------New
            ! CALL DGESV_(m_t,m_t,A,m_t,IPIV,full_Rinv,m_t,INFO )
            ! rones = matmul(B,full_Rinv)
            ! B = rones
            !----------------------------
            If (INFO.ne.0) Then
               Write (6,*) 'k: INFO.ne.0'
               Write (6,*) 'k: INFO=',INFO
               Call Abend()
            End If
            Write (6,*) 'rones ',rones
!
! Now A contains the factors L and U from the factorization A = P*L*U as computed by DGESV
! Where L in the lower triangular matrix with 1 in the diagonal and U is the upper
! triangular matrix of A thus the determinant of A is giving by multipling its diagonal
!
            do i=1,m_t
                diagA(i) = A(i,i)
            enddo
!
!Trend Function (baseline)
            sbO = dot_product(y,B(1:iter))/sum(B(1:iter))
            if (blaAI) then
                sb = y(iter) + blavAI
            else
                if (mblAI) then
                    sb = sbmev
                else
                    if (blAI) Then
                        sb = blvAI
                    else
                        sb = sbO
                    endif
                endif
            endif
            Write (6,*) 'sb ',sb
!
            B = [y-sb,dy]
!
            Ys = B
! ----------------Old calculations --K2
            A=full_r
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            Kv=b
!-----------------New
            ! Kv = matmul(B,full_Rinv)
!------------------------------------
!Likelihood function
            variance = dot_product(Ys,Kv)/m_t
!
            detR = 0
            sign = 1
            do i=1,m_t
                if (diagA(i).le.0) sign=sign*(-1)
                detR = detR + log(abs(diagA(i)))
            enddo
            detR = detR*sign
            lh = variance*exp(detR/dble(m_t))
!
            write(6,*) 'detR',detR
            write(6,*) 'Ys:',Ys
            write(6,*) 'Kv:',Kv
            write(6,*) 'Variance:',variance
            write(6,*) 'm_t',m_t
            write(6,*) 'lh',lh
        END SUBROUTINE k
