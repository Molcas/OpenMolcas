************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine DMINV(N,NMAX,A)
C
C     This subroutine inverts the matrix a of dimension N,N
C
      Implicit real*8 (A-H,O-Z)
      Dimension A(NMAX,NMAX)
      D=1.d0
      Do  K=1,N
        BIGA=A(K,K)
        Do I=1,N
          If (I-K.ne.0) then
           A(I,K)=-A(I,K)/BIGA
          Endif
        Enddo
        Do I=1,N
          HOLD=A(I,K)
          Do J=1,N
           If((I-K)*(J-K).ne.0) then
            A(I,J)=HOLD*A(K,J)+A(I,J)
           Endif
          Enddo
        Enddo
        Do J=1,N
          If(J-K.ne.0) then
           A(K,J)=A(K,J)/BIGA
          Endif
        Enddo
        D=D*BIGA
        A(K,K)=1.D0/BIGA
      Enddo
      Return
      END
