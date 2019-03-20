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
        SUBROUTINE k()
            use globvar
            real*8 B(m_t),A(m_t,m_t) !AF contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
            integer IPIV(m_t),INFO ! ipiv the pivot indices that define the permutation matrix
            B=0
            B(1:iter)=1
            A=full_r !in, coefficent matrix A, out factors L and U from factorization A=PLU on AX=B
            CALL DGESV_(size(A,1), size(shape(B)),A,size(A,2),&
                    IPIV,B,size(B,1),INFO )
            rones=B
            sb=dot_product(y,B(1:iter))/sum(B(1:iter))
            detR=A(1,1)
            do i=2,size(A,1)
                detR=detR*A(i,i)
            enddo
            B=[y-sb,dy]
            Ys=B
            A=full_r
            CALL DGESV_(size(A,1), size(shape(B)),A,size(A,2),&
                    IPIV,B,size(B,1),INFO )
            Kv=b !Kv=K
        !     Write (6,*) 'K',Kv
        END SUBROUTINE k
