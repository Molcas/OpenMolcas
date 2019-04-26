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
            real*8 B(m_t),A(m_t,m_t) !AF contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
            integer IPIV(m_t),INFO,iter ! ipiv the pivot indices that define the permutation matrix
!
!           Initate B according to Eq. (6)
            B=0.0D0
            B(1:iter)=1.0D0
!           Call RecPrt('B',' ',B,1,m_t)
!
!           Initiate A according to Eq. (2)
            A=full_r !in, coefficent matrix A, out factors L and U from factorization A=PLU on AX=B
!           Call RecPrt('A',' ',A,m_t,m_t)
!
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            If (INFO.ne.0) Then
               Write (6,*) 'k: INFO.ne.0'
               Write (6,*) 'k: INFO=',INFO
               Call Abend()
            End If
!
            rones=B
            sb=dot_product(y,B(1:iter))/sum(B(1:iter))
            detR=A(1,1)
            do i=2,m_t
                detR=detR*A(i,i)
            enddo
            B=[y-sb,dy]
            ! write (6,*) 'sb:',sb
            ! write (6,*) 'y',y
            ! write (6,*) 'dy:',dy
            Ys=B
            A=full_r
            CALL DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO )
            Kv=b !Kv=K
!Likelihood function
            ! write(6,*) 'full_R',full_r
            variance=dot_product(Ys,Kv)/m_t
            ! write(6,*) 'detR',detR
            ! write(6,*) 'Ys:',Ys
            ! write(6,*) 'Kv:',Kv
            ! write(6,*) 'Variance:',variance
            lh=variance*exp(log(abs(detr))/m_t)
            ! if (detr<0) then
            !     lh=variance*exp(log(-detr)/m_t)
            ! else
            !     lh=variance*exp(log(detr)/m_t)
            ! endif
        END SUBROUTINE k
