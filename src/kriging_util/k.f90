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
            real*8 B(m_t),A(m_t,m_t),diagA(m_t)
            integer IPIV(m_t),INFO,iter,sign ! ipiv the pivot indices that define the permutation matrix
!
! Initate B according to Eq. (6) of ref.
            B=0.0D0
            B(1:iter)=1.0D0
!           Call RecPrt('B',' ',B,1,m_t)
!
! Initiate A according to Eq. (2) of ref.
            A=full_r !in, coefficent matrix A, out factors L and U from factorization A=PLU on AX=B
!           Call RecPrt('A',' ',A,m_t,m_t)
!
            ! write (6,*) 'A-full_R',A
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            If (INFO.ne.0) Then
               Write (6,*) 'k: INFO.ne.0'
               Write (6,*) 'k: INFO=',INFO
               Call Abend()
            End If
!
! Now A contains the factors L and U from the factorization A = P*L*U as computed by DGESV
! Where L in the lower triangular matrix with 1 in the diagonal and U is the upper
! triangular matrix of A where the determinant is giving by multipling the diagonal
!
            ! sign = 1
            ! detR = 0
            ! do i = 1,m_t
            !     if (A(i,i) < 0) then
            !         sign = sign * (-1)
            !         detR = detR + log(-A(i,i))
            !     else
            !         detR = detR + log(A(i,i))
            !     endif
            !     ! write(6,*) 'detr(i)',i, detr
            !     ! write(6,*) 'a(i,i)',a(i,i)
            ! enddo
            !detR = sign*detR
            do i=1,m_t
                diagA(i) = A(i,i)
            enddo
!
            rones=B
!Trend Function (baseline)
            if (blaAI) then
                sb = y(iter) + blavAI
                sbO = dot_product(y,B(1:iter))/sum(B(1:iter))
            else
                if (mblAI) then
                    sb = sbmev
                    sbO = dot_product(y,B(1:iter))/sum(B(1:iter))
                else
                    if (blAI) Then
                        sb = blvAI
                        sbO = dot_product(y,B(1:iter))/sum(B(1:iter))
                    else
                        sb = dot_product(y,B(1:iter))/sum(B(1:iter))
                    endif
                endif
            endif
!
            B=[y-sb,dy]
            ! write (6,*) 'sb:',sb
            ! write (6,*) 'y',y
            ! write (6,*) 'dy:',dy
            Ys=B
            A=full_r
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            Kv=b !Kv=K
!Likelihood function
            variance = dot_product(Ys,Kv)/m_t
            ! if (maxval(diagA).le.10.0) then
            !     detR = 1!A(1,1)
            !     do i=1,m_t
            !         detR=detR*diagA(i)
            !         write(6,*) 'A, detR:',i, diagA(i), detR
            !     enddo
            !     lh=variance*abs(detR)**(DBLE(1.0D0/DBLE(m_t)))
            ! else
                detR = 0
                sign = 1
                do i=1,m_t
                    if (diagA(i).le.0) sign=sign*(-1)
                    detR = detR + log(abs(diagA(i)))
                    ! write(6,*) 'A, log(detR):',i, diagA(i), detR
                enddo
                ! write(6,*) 'detR/m_t', detR/dble(m_t)
                ! write(6,*) 'exp(detR/m_t)', exp(detR/dble(m_t))
                lh = variance*exp(detR/dble(m_t))
            ! endif
            ! write(6,*) 'detR',detR
            ! write(6,*) 'Ys:',Ys
            ! write(6,*) 'Kv:',Kv
            ! write(6,*) 'Variance:',variance
            ! write(6,*) 'm_t',m_t
            ! lh = variance*exp(detR/m_t)
            !lh=variance*exp(log(detr)/m_t)
            !lh=variance*detR**(DBLE(1.0D0/DBLE(m_t)))
        END SUBROUTINE k
