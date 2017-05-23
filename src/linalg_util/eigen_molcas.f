************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Eigen_Molcas(N,X,D,E)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Compute eigenvalues and eigenvectors of a symmetric real matrix  *
*     by the method of Householder (QR algorithm)                      *
*                                                                      *
*     reference:                                                       *
*     R.S. Martin, C. Reinsch and J.H. Wilkinson                       *
*     Num. Mat. Vol 11, p 181.-195 (1968)                              *
*                                                                      *
*     calling arguments:                                               *
*     N       : Type integer, input.                                   *
*               Dimensions of the matrix X and vectors D and E         *
*     X       : Type real*8 real, input/output                         *
*               on input it is the matrix to diagonalized              *
*               on output it contains the eigenvectors                 *
*     D       : Type real*8 real, output.                              *
*               vector of eigenvalues                                  *
*     E       : Type real*8 real, input/output.                        *
*               Scratch area of length N                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher, University of Lund, Sweden, 2000                 *
*     (the subrotine is based on a old implementation written          *
*      by the comp. center in Munich)                                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension X(N,N),D(N),E(N)

      Parameter ( eps  = 3.0D-17 )
      Parameter ( tol  = 1.0D-22 )
      Parameter ( One  = 1.0D+00 )
      Parameter ( Zero = 0.0D+00 )

      If ( N.eq.1 ) then
        D(1) = X(1,1)
        X(1,1) = One
        Return
      End If
C
C     HOUSEHOLDER REDUCTION
C
      Do ii = 2,N
        i = N+2-ii
        l = i-2
        H = Zero
        G = X(i,i-1)
        If ( l.gt.0 ) then
          Do k = 1,l
            H = H+X(i,k)*X(i,k)
          End Do
          S = H+G*G
          If ( S.lt.Tol ) H = Zero
          If ( H.gt.Zero ) then
            l = l+1
            F = G
            G = sqrt(S)
            If ( G.gt.Zero ) G = -G
            H = S-F*G
            X(i,i-1) = F-G
            F = Zero
            Do j = 1,l
              X(j,i) = X(i,j)/H
              sum = Zero
              Do k = 1,j
                sum = sum+X(j,k)*X(i,k)
              End Do
              Do k = j+1,l
                sum = sum+X(k,j)*X(i,k)
              End Do
              E(j) = sum/H
              F = F+sum*X(j,i)
            End Do
            F = F/(H+H)
            Do j = 1,l
              E(j) = E(j)-F*X(i,j)
            End Do
            Do j = 1,l
              F = X(i,j)
              scal = E(j)
              Do k = 1,j
                X(j,k) = X(j,k)-F*E(k)-X(i,k)*scal
              End Do
            End Do
          End If
        End If
        D(i) = H
        E(i-1) = G
      End Do
C
C     ACCUMULATION OF TRANSFORMATION MATRICES
C
      D(1) = X(1,1)
      X(1,1) = One
      Do i = 2,N
        l = i-1
        If ( D(i).gt.0.0D0 ) then
          Do j = 1,l
            S = Zero
            Do k = 1,l
              S = S+X(i,k)*X(k,j)
            End Do
            Do k = 1,l
              X(k,j) = X(k,j)-S*X(k,i)
            End Do
          End Do
        End If
        D(i) = X(i,i)
        X(i,i) = One
        Do j = 1,l
          X(i,j) = Zero
          X(j,i) = Zero
        End Do
      End Do
C
C     DIAGONAlIZATION OF THE TRIDIAGONAL MATRIX
C
      B = Zero
      F = Zero
      E(N) = Zero
      Do l = 1,N
        H = eps*(abs(D(l))+abs(E(l)))
        If ( H.GT.B ) B = H
        Do j = l,N
          If ( abs(E(j)).lE.B ) then
            If ( j.eq.l ) Goto 100
          End If
        End Do
        j = N
        Do while ( abs(E(l)).GT.B )
          P = (D(l+1)-D(l))*0.5D0/E(l)
          R = sqrt(P*P+One)
          If ( P.ge.Zero ) then
            P = P+R
          Else
            P = P-R
          End If
          H = D(l)-E(l)/P
          Do i = l,N
            D(i) = D(i)-H
          End Do
          F = F+H
          P = D(j)
          C = One
          S = Zero
          Do ii = l,j-1
            i = l+j-1-ii
            G = C*E(i)
            H = C*P
            If ( abs(P).lT.abs(E(i)) ) then
              C = P/E(i)
              R = sqrt(C*C+One)
              E(i+1) = S*E(i)*R
              S = One/R
              C = C/R
            Else
              C = E(i)/P
              R = sqrt(C*C+One)
              E(i+1) = S*P*R
              S = C/R
              C = One/R
            End If
            P = C*D(i)-S*G
            D(i+1) = H+S*(C*G+S*D(i))
            Do k = 1,N
              H = X(k,i+1)
              X(k,i+1) = X(k,i)*S+H*C
              X(k,i) = X(k,i)*C-H*S
            End Do
          End Do
          E(l) = S*P
          D(l) = C*P
        End Do
100     Continue
        D(l) = D(l)+F
      End Do
C
C     ORDERING OF EIGENVALUES
C
      Do i = 1,N-1
        Do j = i+1,N
          If ( D(j).lt.D(i) ) then
            tmp = D(j)
            D(j) = D(i)
            D(i) = tmp
            Do k = 1,N
              tmp = X(k,j)
              X(k,j) = X(k,i)
              X(k,i) = tmp
            End Do
          End If
        End Do
      End Do
C
C     FIXING OF SIGN
C
      Do i = 1,N
        k = 1
        tmp = abs(X(k,i))
        Do j = 2,N
          If ( tmp.le.abs(X(j,i)) ) then
            tmp = abs(X(j,i))
            k = j
          End If
        End Do
        If ( X(k,i).lt.Zero ) then
          Do j = 1,N
            X(j,i) = -X(j,i)
          End Do
        End If
      End Do

      Return
      End
