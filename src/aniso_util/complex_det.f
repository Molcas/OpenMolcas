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
      Complex*16 FUNCTION CDET(MA,N,A)
C================================================= ========
c MA is the MAximal dimension
c N is the dimension N<MA
c A is the Complex matrix
c================================================= ======
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)            :: N,MA
      Complex(kind=8), intent(inout):: A(MA,MA)
      Integer                        :: I,J,K,L,K1
      Real(kind=8)                  :: P,Q
      Complex(kind=8)               :: CP,CQ
C
      I=0
      CDET = (0.0_wp,0.0_wp)

      Do K=1,N
      P = 0.0_wp

         Do J=K,N
            Q = ABS(A(J,K))
            If (Q.gt.P) Then
               P = Q
               I = J
            End If
         End Do

         CP = (1.0_wp,0.0_wp)/A(I,K)

         If (I.ne.K) Then
            CDET = -CDET
            Do L=K,N
               CQ = A(I,L)
               A(I,L) = A(K,L)
               A(K,L) = CQ
            End Do
         End If

         CDET = CDET*A(K,K)

         If (K.lt.N) Then
            K1 = K+1
            Do I = K1,N
               CQ=A(I,K)*CP
               Do L = K1,N
                  A(I,L) = A(I,L)-CQ*A(K,L)
               End Do
            End Do
         End If

      End Do ! K

      Return
      End
