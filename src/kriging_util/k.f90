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
#include "stdalloc.fh"
            real*8, Allocatable:: B(:), A(:,:)
            Integer, Allocatable:: IPIV(:)
            integer INFO,iter,sign ! ipiv the pivot indices that define the permutation matrix
!
            Call mma_Allocate(B,m_t,Label="B")
            Call mma_Allocate(A,m_t,m_t,Label="A")
            Call mma_Allocate(IPIV,m_t,Label="IPIV")
!
! Initiate B according to Eq. (6) of ref.
            B=0.0D0
            B(1:iter)=1.0D0
!
! Initiate A according to Eq. (2) of ref.
!
            !-----------------New
            ! Full_RInv=0
            ! forall(j=1:m_t) Full_RInv(j,j)=1
            ! A = full_r !in, coefficent matrix A, out factors L and U from factorization A=PLU on AX=B
            ! CALL DGESV_(m_t,m_t,A,m_t,IPIV,full_Rinv,m_t,INFO )
            ! rones = matmul(full_Rinv,B)!matmul(B,full_Rinv)
            ! B = rones
            ! Write (6,*) 'A new ',A
            ! ----------------Old calculations --K1
            A(:,:) = full_r(:,:)
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            rones(:)=B(:)
            ! Write (6,*) 'A old ',A
            !----------------------------
            If (INFO.ne.0) Then
               Write (6,*) 'k: INFO.ne.0'
               Write (6,*) 'k: INFO=',INFO
               Call Abend()
            End If
!
! Now A contains the factors L and U from the factorization A = P*L*U as computed by DGESV
! Where L in the lower triangular matrix with 1 in the diagonal and U is the upper
! triangular matrix of A thus the determinant of A is giving by multipling its diagonal
!
            detR = 0.0d0
            sign = 1
            do i=1,m_t
                ! if (A(i,i).le.0) sign=sign*(-1)
                detR = detR + log(abs(A(i,i)))
            enddo
            ! detR = detR*sign
!
!Trend Function (baseline)
            sbO = dot_product(y,B(1:iter))/sum(B(1:iter))
            if (blaAI) then
                sb = y(iter) + blavAI
            else if (mblAI) then
                sb = sbmev
            else if (blAI) Then
                sb = blvAI
            else
                 sb = sbO
            endif
        !   Write (6,*) 'K: sb=',sb
        !   Write (6,*) 'K: y=',y
        !   Write (6,*) 'K: dy=',dy
!
            B(:) = [y-sb,dy]
            Kv(:)=B
! ----------------Old calculations --K2
            A(:,:)=full_r
!           Write (6,*) 'K: y=',y
        !   Write (6,*) 'K: B(Ys)=',B
!           Write (6,*) 'K: A=',A
            CALL DGESV_(m_t,1,A,m_t,IPIV,Kv,m_t,INFO)
!-----------------New
            ! Kv = matmul(B,full_Rinv)
!------------------------------------
!Likelihood function
            variance = dot_product(B,Kv)/m_t
            lh = variance*exp(detR/dble(m_t))
!
            ! write(6,*) 'detR',detR
            ! write(6,*) 'Kv orig:',Kv
            ! write(6,*) 'Variance:',variance
            ! write(6,*) 'm_t',m_t
            ! write(6,*) 'lh',lh
!
            Call mma_Deallocate(B)
            Call mma_Deallocate(A)
            Call mma_Deallocate(IPIV)
        END SUBROUTINE k