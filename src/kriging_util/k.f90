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
            real*8, Allocatable:: B(:), A(:,:), diagA(:), iden(:,:), Ys(:)
            Integer, Allocatable:: IPIV(:)
            integer INFO,iter,sign ! ipiv the pivot indices that define the permutation matrix
!
            Call mma_Allocate(B,m_t,Label="B")
            Call mma_Allocate(A,m_t,m_t,Label="A")
            Call mma_Allocate(diagA,m_t,Label="diagAB")
            Call mma_Allocate(iden,m_t,m_t,Label="iden")
            Call mma_Allocate(IPIV,m_t,Label="IPIV")
            Call mma_Allocate(Ys,m_t,Label="Ys")
!
! Initiate B according to Eq. (6) of ref.
            B=0.0D0
            B(1:iter)=1.0D0
!
            call miden(iden,m_t)
            full_Rinv = iden
! Initiate A according to Eq. (2) of ref.
!
            !-----------------New
            ! A = full_r !in, coefficent matrix A, out factors L and U from factorization A=PLU on AX=B
            ! CALL DGESV_(m_t,m_t,A,m_t,IPIV,full_Rinv,m_t,INFO )
            ! rones = matmul(full_Rinv,B)!matmul(B,full_Rinv)
            ! B = rones
            ! Write (6,*) 'A new ',A
            ! ----------------Old calculations --K1
            A = full_r
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            rones=B
            ! Write (6,*) 'A old ',A
            !----------------------------
            If (INFO.ne.0) Then
               Write (6,*) 'k: INFO.ne.0'
               Write (6,*) 'k: INFO=',INFO
               Call Abend()
            End If
            ! Call RecPrt('full_Rinv',  ' ',full_Rinv,m_t,m_t)
            ! Call RecPrt('full_Rinv*full_R = I',  ' ',matmul(full_Rinv,full_r),m_t,m_t)
            ! Write (6,*) 'rones ',rones
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
            else if (mblAI) then
                sb = sbmev
            else if (blAI) Then
                sb = blvAI
            else
                 sb = sbO
            endif
!           Write (6,*) 'K: sb=',sb
!           Write (6,*) 'K: y=',y
!           Write (6,*) 'K: dy=',dy
!
            B = [y-sb,dy]
!
            Ys = B
! ----------------Old calculations --K2
            A=full_r
!           Write (6,*) 'K: y=',y
!           Write (6,*) 'K: B=',B
!           Write (6,*) 'K: A=',A
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO)
            Kv=b
!-----------------New
            ! Kv = matmul(B,full_Rinv)
!------------------------------------
!Likelihood function
            variance = dot_product(Ys,Kv)/m_t
!
            detR = 0.0d0
            sign = 1
            do i=1,m_t
                ! if (diagA(i).le.0) sign=sign*(-1)
                detR = detR + log(abs(diagA(i)))
            enddo
            ! detR = detR*sign
            lh = variance*exp(detR/dble(m_t))
!
            ! write(6,*) 'detR',detR
            ! write(6,*) 'Ys:',Ys
            ! write(6,*) 'Kv:',Kv
            ! write(6,*) 'Variance:',variance
            ! write(6,*) 'm_t',m_t
            ! write(6,*) 'lh',lh
!
            Call mma_Deallocate(B)
            Call mma_Deallocate(A)
            Call mma_Deallocate(diagA)
            Call mma_Deallocate(iden)
            Call mma_Deallocate(IPIV)
            Call mma_Deallocate(Ys)
        END SUBROUTINE k
