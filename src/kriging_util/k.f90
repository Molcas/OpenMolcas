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
!
            B(:) = [y-sb,dy]
            Kv(:)=B
!
            A(:,:)=full_r
!
            CALL DGESV_(m_t,1,A,m_t,IPIV,Kv,m_t,INFO)
!
!Likelihood function
            variance = dot_product(B,Kv)/m_t
            lh = variance*exp(detR/dble(m_t))
!
            Call mma_Deallocate(B)
            Call mma_Deallocate(A)
            Call mma_Deallocate(IPIV)
        END SUBROUTINE k